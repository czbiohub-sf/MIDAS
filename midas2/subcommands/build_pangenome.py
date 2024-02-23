#!/usr/bin/env python3
import os
import sys
import json
from collections import defaultdict
from multiprocessing import Semaphore
import Bio.SeqIO
from midas2.common.argparser import add_subcommand
from midas2.common.utils import tsprint, InputStream, OutputStream, retry, command, hashmap, split, multiprocessing_map, multithreading_map, num_vcpu, select_from_tsv, transpose, find_files, upload, upload_star, flatten, pythonpath, num_physical_cores, copy_star
from midas2.common.utilities import decode_species_arg, has_ambiguous_bases
from midas2.models.midasdb import MIDAS_DB
from midas2.params.inputs import MIDASDB_NAMES


"""
Input:
    - Cleaned Prokka annotations of all the genomes for given species
Output:
    - temp/vsearch/genes.ffn
    - temp/vsearch/genes.len
    - temp/vsearch/gene_info.txt
    - temp/vsearch/centroid.xx.ffn # with --recluster, we only have xx = 99
"""


CLUSTERING_PERCENTS = [99, 95, 90, 85, 80, 75]
CLUSTERING_PERCENTS = sorted(CLUSTERING_PERCENTS, reverse=True)


# Up to this many concurrent species builds.
CONCURRENT_SPECIES_BUILDS = Semaphore(3)


def destpath(midas_db, species_id, filename):
    return midas_db.get_target_layout("pangenome_file", True, species_id, "", filename)

def localpath(midas_db, species_id, filename):
    return midas_db.get_target_layout("pangenome_file", False, species_id, "", filename)

def localtemp(midas_db, species_id, step, filename):
    return midas_db.get_target_layout("pangenome_tempfile", False, species_id, step, filename)


@retry
def find_files_with_retry(f):
    return find_files(f)


# 1. Occasional failures in aws s3 cp require a retry.
# 2. In future, for really large numbers of genomes, we may prefer a separate wave of retries for all first-attempt failures.
# 3. The Bio.SeqIO.parse() code is CPU-bound and thus it's best to run this function in a separate process for every genome.
@retry
def clean_genes(packed_ids):
    """
    1. extra newlines have been removed from the gene sequences, so that each gene sequence occupies a single line below the gene header
    2. DNA letters have been converted to uppercase, and
    3. degenerate genes (with empty sequences or headers containing "|") have been excluded
    4. Remove genes with ambiguous bases
    5. remove short genes (<= 200 bps)
    """
    genome_id, genes_ffn = packed_ids

    output_genes = f"{genome_id}.genes.ffn"
    output_len = f"{genome_id}.genes.len"
    with open(output_genes, 'w') as o_genes, \
         open(output_len, 'w') as o_info, \
         InputStream(genes_ffn, check_path=False) as genes:  # check_path=False because for flat directory structure it's slow
        for rec in Bio.SeqIO.parse(genes, 'fasta'):
            gene_id = rec.id
            gene_seq = str(rec.seq).upper()
            gene_len = len(gene_seq)
            if gene_len <= 200 or has_ambiguous_bases(gene_seq) or gene_id == '' or gene_id == '|':
                pass
            else:
                o_genes.write(f">{gene_id}\n{gene_seq}\n")
                o_info.write(f"{gene_id}\t{genome_id}\t{gene_len}\n")
    return output_genes, output_len


def clean_centroids(percent_id, centroids_fp):
    # Check ambiguous centroids.99, and write to separate files if exist.
    output_ambiguous = f"centroids.{percent_id}.ambiguous.ffn"
    output_clean = f"centroids.{percent_id}.clean.ffn"
    with open(output_ambiguous, 'w') as o_ambiguous, \
         open(output_clean, 'w') as o_clean, \
         InputStream(centroids_fp, check_path=False) as centroids:  # check_path=False because for flat directory structure it's slow
        for rec in Bio.SeqIO.parse(centroids, 'fasta'):
            c_id = rec.id
            c_seq = str(rec.seq).upper()
            if has_ambiguous_bases(c_seq):
                o_ambiguous.write(f">{c_id}\n{c_seq}\n")
            else:
                o_clean.write(f">{c_id}\n{c_seq}\n")
    return output_ambiguous, output_clean


def vsearch(percent_id, genes, num_threads=num_vcpu):
    centroids = f"centroids.{percent_id}.ffn"
    uclust = f"uclust.{percent_id}.txt"
    if find_files(centroids) and find_files(uclust):
        tsprint(f"Found vsearch results at percent identity {percent_id} from prior run.")
    else:
        try:
            command(f"vsearch --quiet --cluster_fast {genes} --id {percent_id/100.0} --threads {num_threads} --centroids {centroids} --uc {uclust}")
        except:
            # Do not keep bogus zero-length files;  those are harmful if we rerun in place.
            command(f"mv {centroids} {centroids}.bogus", check=False)
            command(f"mv {uclust} {uclust}.bogus", check=False)
            raise
    return centroids, uclust


def parse_uclust(uclust_file, select_columns):
    # The uclust TSV file does not contain a header line.  So, we have to hardcode the schema here.  Then select specified columns.
    all_uclust_columns = ['type', 'cluster_id', 'size', 'pid', 'strand', 'skip1', 'skip2', 'skip3', 'gene_id', 'centroid_id']
    with InputStream(uclust_file) as ucf:
        for r in select_from_tsv(ucf, select_columns, all_uclust_columns):
            yield r


def read_uclust_info(centroid_info, uclust_file, percent_id):
    # Get centroid_info from uclust
    for r_type, r_gene, r_centroid in parse_uclust(uclust_file, ['type', 'gene_id', 'centroid_id']):
        if r_type == 'S':
            # r itself is the centroid of its cluster
            centroid_info[r_gene][percent_id] = r_gene
        elif r_type == 'H':
            # r is not itself a centroid
            centroid_info[r_gene][percent_id] = r_centroid
        else:
            # ignore all other r types
            pass


def xref(cluster_files):
    """
    Produce the per-species gene_info.txt file
    """
    # Let centroid_info[gene][percent_id] be the centroid of the percent_id cluster contianing gene.
    # The max_percent_id centroids are computed directly for all genes.  Only these centroids are
    # then reclustered to lower percent_id's.
    #
    # The centroids are themselves genes, and their ids, as all gene_ids, are strings
    # generated by the annotation tool Prokka.
    centroid_info = defaultdict(dict)
    for percent_id, (_, uclust_file) in cluster_files.items():
        read_uclust_info(centroid_info, uclust_file, percent_id)

    # Check for a problem that occurs with improper import of genomes (when contig names clash).
    percents = cluster_files.keys()
    max_percent_id = max(percents)
    for g in centroid_info:
        cg = centroid_info[g][max_percent_id]
        ccg = centroid_info[cg][max_percent_id]
        assert cg == ccg, f"The {max_percent_id}-centroid relation should be idempotent, however, {cg} != {ccg}.  See https://github.com/czbiohub/MIDAS2/issues/16"

    # At this point we have the max_percent_id centroid for any gene gc, but we lack
    # coarser clustering assignments for many genes -- we only have those for genes
    # that are themelves centroids of max_percent_id clusters.
    #
    # We can infer the remaining cluster assignments for all genes by transitivity.
    # For any gene gc, look up the clusters containing gc's innermost centroid,
    # gc[max_percent_id].  Those clusters also contain gc.
    for gc in centroid_info.values():
        gc_recluster = centroid_info[gc[max_percent_id]]
        for percent_id in percents:
            gc[percent_id] = gc_recluster[percent_id]

    return centroid_info


def write_gene_info(centroid_info, percents, gene_info_file):
    # Write centroid_info[gene][percent_id] to gene_info_file
    with OutputStream(gene_info_file) as gene_info:
        header = ['gene_id'] + [f"centroid_{pid}" for pid in percents]
        gene_info.write('\t'.join(header) + '\n')
        genes = centroid_info.keys()
        for gene_id in sorted(genes):
            gene_info.write(gene_id)
            for centroid in centroid_info[gene_id].values():
                gene_info.write('\t')
                gene_info.write(centroid)
            gene_info.write('\n')


def build_pangenome(args):
    if args.zzz_worker_mode:
        build_pangenome_worker(args)
    else:
        build_pangenome_master(args)


def build_pangenome_master(args):

    # Fetch table of contents from s3.
    # This will be read separately by each species build subcommand, so we make a local copy.
    midas_db = MIDAS_DB(os.path.abspath(args.midasdb_dir), args.midasdb_name)
    species = midas_db.uhgg.species
    num_threads = args.num_threads

    def species_work(species_id):
        assert species_id in species, f"Species {species_id} is not in the database."
        species_genomes = species[species_id]

        # The species build will upload this file last, after everything else is successfully uploaded.
        # Therefore, if this file exists in s3, there is no need to redo the species build.
        dest_file = destpath(midas_db, species_id, "temp/vsearch/gene_info.txt") # TODO: is this still our final dest file?
        local_file = localpath(midas_db, species_id, "temp/vsearch/gene_info.txt")

        msg = f"Building pangenome for species {species_id} with {len(species_genomes)} total genomes."
        if args.upload and find_files_with_retry(dest_file):
            if not args.force:
                tsprint(f"Destination {dest_file} for species {species_id} pangenome already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Building", "Rebuilding")
        if not args.upload and os.path.exists(local_file):
            if not args.force:
                tsprint(f"Destination {local_file} for species {species_id} pangenome already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Building", "Rebuilding")

        with CONCURRENT_SPECIES_BUILDS:
            tsprint(msg)
            last_dest = midas_db.get_target_layout("pangenome_log", True, species_id)
            local_dest = midas_db.get_target_layout("pangenome_log", False, species_id)
            local_dir = os.path.dirname(local_dest)
            command(f"mkdir -p {local_dir}/temp/vsearch")

            worker_log = os.path.basename(local_dest)
            worker_subdir = os.path.dirname(local_dest) if args.scratch_dir == "." else f"{args.scratch_dir}/buildpan/{species_id}"
            worker_log = f"{worker_subdir}/{worker_log}"
            if not args.debug:
                command(f"rm -rf {worker_subdir}")

            worker_subdir = f"{worker_subdir}/temp/vsearch"
            if not os.path.isdir(worker_subdir):
                command(f"mkdir -p {worker_subdir}")

            # Recurisve call via subcommand.  Use subdir, redirect logs.
            subcmd_str = f"--zzz_worker_mode -t {num_threads} --midasdb_name {args.midasdb_name} --midasdb_dir {os.path.abspath(args.midasdb_dir)} {'--debug' if args.debug else ''} {'--upload' if args.upload else ''} --scratch_dir {args.scratch_dir}"
            worker_cmd = f"cd {worker_subdir}; PYTHONPATH={pythonpath()} {sys.executable} -m midas2 build_pangenome -s {species_id} {subcmd_str} {'--recluster' if args.recluster else ''} &>> {worker_log}"
            with open(f"{worker_log}", "w") as slog:
                slog.write(msg + "\n")
                slog.write(worker_cmd + "\n")

            try:
                command(worker_cmd)
            finally:
                # Cleanup should not raise exceptions of its own, so as not to interfere with any
                # prior exceptions that may be more informative.  Hence check=False.
                if args.upload:
                    upload(f"{worker_log}", last_dest, check=False)
                if args.scratch_dir != ".":
                    command(f"cp -r {worker_log} {local_dest}")
                # Clean up temporary files
                if not args.debug:
                    command(f"rm -rf {worker_subdir}", check=False)

    # Check for destination presence in s3 with up to 10-way concurrency.
    # If destination is absent, commence build with up to 3-way concurrency as constrained by CONCURRENT_SPECIES_BUILDS.
    species_id_list = decode_species_arg(args, species)
    multithreading_map(species_work, species_id_list, num_threads=num_threads)


def build_pangenome_worker(args):

    violation = "Please do not call build_pangenome_worker directly.  Violation"
    assert args.zzz_worker_mode, f"{violation}:  Missing --zzz_worker_mode arg."

    species_id = args.species
    midas_db = MIDAS_DB(args.midasdb_dir, args.midasdb_name, args.num_threads) #<----

    species = midas_db.uhgg.species
    assert species_id in species, f"build_pangenome_worker::Species {species_id} is not in the database."

    ffns_by_genomes = midas_db.fetch_files("annotation_ffn", [species_id], False)
    cleaned = multiprocessing_map(clean_genes, ((genome_id, geneffn) for genome_id, geneffn in ffns_by_genomes.items()), args.num_threads)

    command("rm -f genes.ffn genes.len")
    for temp_files in split(cleaned, 20):  # keep "cat" commands short
        ffn_files, len_files = transpose(temp_files)
        command("cat " + " ".join(ffn_files) + " >> genes.ffn")
        command("cat " + " ".join(len_files) + " >> genes.len")
        command("rm " + " ".join(ffn_files))
        command("rm " + " ".join(len_files))

    # The initial clustering to max_percent takes longest.
    max_percent, lower_percents = CLUSTERING_PERCENTS[0], CLUSTERING_PERCENTS[1:]
    cluster_files = {max_percent: vsearch(max_percent, "genes.ffn")}

    # Check ambiguous bases of centroids.99
    c99_ambigous, c99_clean = clean_centroids(max_percent, cluster_files[max_percent][0])

    # TODO: we only implement recluster on HPC; and we highly recommend recluster
    if not args.recluster:
        # Reclustering of the max_percent centroids is usually quick, and can proceed in prallel.
        recluster = lambda percent_id: vsearch(percent_id, cluster_files[max_percent][0])
        cluster_files.update(hashmap(recluster, lower_percents))

    centroid_info = xref(cluster_files)
    write_gene_info(centroid_info, cluster_files.keys(), "gene_info.txt")

    if args.upload:
        # Create list of (source, dest) pairs for uploading.
        # Note that centroids.{max_percent}.ffn is uploaded to two different destinations.
        upload_tasks = [
            ("genes.ffn", destpath(midas_db, species_id, "temp/vsearch/genes.ffn")),
            ("genes.len", destpath(midas_db, species_id, "temp/vsearch/genes.len")),
        ]
        for src in flatten(cluster_files.values()):
            upload_tasks.append((src, destpath(midas_db, species_id, f"temp/vsearch/{src}"))) #<---
        command(f"aws s3 rm --recursive {os.path.dirname(last_dest_file)}")
        # Upload in parallel.
        multithreading_map(upload_star, upload_tasks)
        # Leave this upload for last, so the presence of this file in s3 would indicate the entire species build has succeeded.
        last_output = "gene_info.txt"
        last_dest_file = destpath(midas_db, species_id, f"temp/vsearch/{last_output}")
        upload(last_output, last_dest_file)

    if args.scratch_dir != ".":
        copy_tasks = [
            ("genes.ffn", localtemp(midas_db, species_id, "vsearch", "genes.ffn")),
            ("genes.len", localtemp(midas_db, species_id, "vsearch", "genes.len")),
            ("gene_info.txt", localtemp(midas_db, species_id, "vsearch", "gene_info.txt")),
            (c99_ambigous, localtemp(midas_db, species_id, "vsearch", c99_ambigous)),
            (c99_clean, localtemp(midas_db, species_id, "vsearch", c99_clean))
        ]
        for src in flatten(cluster_files.values()):
            copy_tasks.append((src, localtemp(midas_db, species_id, "vsearch", src)))
        multithreading_map(copy_star, copy_tasks, args.num_threads)


def register_args(main_func):
    subparser = add_subcommand('build_pangenome', main_func, help='Build pangenome for specified species')
    subparser.add_argument('-s',
                           '--species',
                           dest='species',
                           required=False,
                           help="species[,species...] whose pangenome(s) to build;  alternatively, species slice in format idx:modulus, e.g. 1:30, meaning build species whose ids are 1 mod 30; or, the special keyword 'all' meaning all species")
    subparser.add_argument('--midasdb_name',
                           dest='midasdb_name',
                           type=str,
                           default="uhgg",
                           choices=MIDASDB_NAMES,
                           help="MIDAS Database name.")
    subparser.add_argument('--midasdb_dir',
                           dest='midasdb_dir',
                           type=str,
                           default=".",
                           help="Path to local MIDAS Database.")
    subparser.add_argument('-t',
                           '--num_threads',
                           dest='num_threads',
                           type=int,
                           default=num_physical_cores,
                           help="Number of threads")
    subparser.add_argument('--upload',
                           action='store_true',
                           default=False,
                           help="Upload built files to AWS S3.")
    subparser.add_argument('--scratch_dir',
                           dest='scratch_dir',
                           type=str,
                           default=".",
                           help="Absolute path to scratch directory for fast I/O.")
    subparser.add_argument('--recluster',
                           action='store_true',
                           default=False,
                           help="Recluster centroids.99 via cd-hit.")
    return main_func


@register_args
def main(args):
    tsprint(f"Build species-specific pangenome in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    build_pangenome(args)
