#!/usr/bin/env python3
import os
import sys
from collections import defaultdict
from multiprocessing import Semaphore

from midas2.common.argparser import add_subcommand
from midas2.common.utils import tsprint, retry, InputStream, OutputStream, command, hashmap, cat_files, multithreading_map, select_from_tsv, pythonpath, num_physical_cores, copy_star, flatten
from midas2.common.utilities import decode_species_arg
from midas2.models.midasdb import MIDAS_DB
from midas2.params.inputs import MIDASDB_NAMES
from midas2.params.schemas import GENE_LENGTH_SCHEMA, MARKER_INFO_SCHEMA, PANGENOME_INFO_SCHEMA
from midas2.subcommands.build_pangenome import vsearch, CLUSTERING_PERCENTS, read_uclust_info, localpath


"""
Input:
    - temp/cdhit/gene_info.txt
    - temp/cdhit/centroids.99.ffn
    - temp/cdhit/genes.ffn
    - temp/cdhit/genes.len
Output:
    - genes_info.tsv: the unique id of this table is gene id of the cleaned set of genes
    - centroids.99.ffn
    - temp/centroids.xx.ffn
"""


# Up to this many concurrent species builds.
CONCURRENT_SPECIES_BUILDS = Semaphore(1)


@retry
def scan_gene_length(gene_length_file):
    gene_length_dict = {}
    with InputStream(gene_length_file) as stream:
        for r in select_from_tsv(stream, schema=GENE_LENGTH_SCHEMA, result_structure=dict):
            gene_length_dict[r["gene_id"]] = r["gene_length"]
    return gene_length_dict


@retry
def scan_mapfile(mapfile):
    """ Extract <marker_id, gene_id> pairs from the marker mapfile """
    dict_of_markers = {}
    with InputStream(mapfile) as stream:
        # Example: 100001	GUT_GENOME000001	GUT_GENOME000001_01635	906	B000032
        for gene_id, marker_id in select_from_tsv(stream, ["gene_id", "marker_id"], schema=MARKER_INFO_SCHEMA):
            dict_of_markers[gene_id] = marker_id
    return dict_of_markers


def read_gene_info(centroid_info, gene_info_file, percent_id):
    # Parse intermediate gene_info_cdhit.tsv
    with InputStream(gene_info_file) as stream:
        for r in select_from_tsv(stream, selected_columns=['gene_id', 'centroid_99'], schema={'gene_id':str, 'centroid_99':str}):
            centroid_info[r[0]][percent_id] = r[1]


def augment_gene_info(centroid_info, gene_to_marker, dict_of_gene_length, gene_info_file):
    """ Augment gene_info.txt with two additional columns: gene_length and marker_id """
    with OutputStream(gene_info_file) as stream:
        stream.write("\t".join(PANGENOME_INFO_SCHEMA.keys()) + "\n")
        for gene_id, r in centroid_info.items():
            gene_len = dict_of_gene_length[gene_id]
            marker_id = gene_to_marker[gene_id] if gene_id in gene_to_marker else ""
            val = [gene_id] + list(r.values()) + [gene_len, marker_id]
            stream.write("\t".join(map(str, val)) + "\n")


def xref(cluster_files):
    # Again, let centroid_info[gene][percent_id] be the centroid of the percent_id
    # cluster containing gene; then reclustered to lower percent_id's.
    # Output: gene_info.txt for filtered centroids_99 genereated by cd-hit

    percents = cluster_files.keys()
    max_percent_id = max(percents)
    centroid_info = defaultdict(dict)
    for percent_id, (_, uclust_file) in cluster_files.items():
        if percent_id == max_percent_id:
            read_gene_info(centroid_info, uclust_file, percent_id)
        else:
            read_uclust_info(centroid_info, uclust_file, percent_id)

    # Check for a problem that occurs with improper import of genomes (when contig names clash).
    for g in centroid_info:
        cg = centroid_info[g][max_percent_id]
        ccg = centroid_info[cg][max_percent_id]
        assert cg == ccg, f"The {max_percent_id}-centroid relation should be idempotent, however, {cg} != {ccg}."

    # Infer coarser clusters assignments for all genes by transitivity
    for gc in centroid_info.values():
        gc_recluster = centroid_info[gc[max_percent_id]]
        for percent_id in percents:
            gc[percent_id] = gc_recluster[percent_id]

    return centroid_info


def recluster_centroid(args):
    if args.zzz_worker_mode:
        recluster_centroid_worker(args)
    else:
        recluster_centroid_master(args)


def recluster_centroid_master(args):

    # Fetch table of contents from s3.
    # This will be read separately by each species build subcommand, so we make a local copy.
    midas_db = MIDAS_DB(os.path.abspath(args.midasdb_dir), args.midasdb_name)
    species = midas_db.uhgg.species
    num_threads = args.num_threads

    def species_work(species_id):
        assert species_id in species, f"Species {species_id} is not in the database."
        species_genomes = species[species_id]

        # The species build will upload this file last, after everything else is successfully uploaded.
        # Therefore, if this file exists locally, there is no need to redo the species build.
        dest_file = midas_db.get_target_layout("pangenome_genes_info", False, species_id)
        msg = f"Reclustering centroids for species {species_id} with {len(species_genomes)} total genomes."

        if os.path.exists(dest_file):
            if not args.force:
                tsprint(f"Destination {dest_file} for species {species_id} pangenome already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Building", "Rebuilding")

        with CONCURRENT_SPECIES_BUILDS:
            tsprint(msg)
            last_dest = midas_db.get_target_layout("recluster_log", True, species_id)
            local_dest = midas_db.get_target_layout("recluster_log", False, species_id)
            local_dir = os.path.dirname(local_dest)
            command(f"mkdir -p {local_dir}/temp/cdhit")

            worker_log = os.path.basename(local_dest)
            worker_subdir = os.path.dirname(local_dest) if args.scratch_dir == "." else f"{args.scratch_dir}/buildpan/{species_id}"
            worker_log = f"{worker_subdir}/{worker_log}"
            if not args.debug:
                command(f"rm -rf {worker_subdir}")

            worker_subdir = f"{worker_subdir}/temp"
            if not os.path.isdir(f"{worker_subdir}/cdhit"):
                command(f"mkdir -p {worker_subdir}/cdhit")

            # Recurisve call via subcommand.  Use subdir, redirect logs.
            subcmd_str = f"--zzz_worker_mode -t {num_threads} --midasdb_name {args.midasdb_name} --midasdb_dir {os.path.abspath(args.midasdb_dir)} {'--debug' if args.debug else ''} {'--upload' if args.upload else ''} --scratch_dir {args.scratch_dir}"
            worker_cmd = f"cd {worker_subdir}; PYTHONPATH={pythonpath()} {sys.executable} -m midas2 recluster_centroids -s {species_id} {subcmd_str} &>> {worker_log}"

            with open(f"{worker_log}", "w") as slog:
                slog.write(msg + "\n")
                slog.write(worker_cmd + "\n")

            try:
                command(worker_cmd)
            finally:
                if not args.debug:
                    command(f"rm -rf {worker_subdir}", check=False)
                if args.scratch_dir != ".":
                    command(f"cp -r {worker_log} {local_dest}")
                # TODO: handle upload

    # Check for destination presence in s3 with up to 8-way concurrency.
    # If destination is absent, commence build with up to 3-way concurrency as constrained by CONCURRENT_SPECIES_BUILDS.
    species_id_list = decode_species_arg(args, species)
    CONCURRENT_RUNS = min(1, int(args.num_threads / 8))
    multithreading_map(species_work, species_id_list, num_threads=CONCURRENT_RUNS)


def recluster_centroid_worker(args):

    violation = "Please do not call generate_gene_info_worker directly.  Violation"
    assert args.zzz_worker_mode, f"{violation}:  Missing --zzz_worker_mode arg."

    species_id = args.species
    midas_db = MIDAS_DB(args.midasdb_dir, args.midasdb_name)

    species = midas_db.uhgg.species
    assert species_id in species, f"{violation}: Species {species_id} is not in the database."

    # start with existing centroids_99 generated by cd-hit
    max_percent, lower_percents = CLUSTERING_PERCENTS[0], CLUSTERING_PERCENTS[1:]
    cluster_files = {max_percent: ("cdhit/centroids.99.ffn", "cdhit/gene_info.txt")}
    assert os.path.isfile("cdhit/centroids.99.ffn"), "Error: need to run cd-hit denoise pipeline."
    assert os.path.isfile("cdhit/gene_info.txt"), "Error: need to run cd-hit denoise pipeline."

    # reclustering of the max_percent centroids is usually quick, and can proceed in prallel.
    recluster = lambda percent_id: vsearch(percent_id, cluster_files[max_percent][0], num_threads=args.num_threads)
    cluster_files.update(hashmap(recluster, lower_percents))

    # cluster centroids.99 to lower level
    centroid_info = xref(cluster_files)

    # augment temp/vsearch/gene_info.txt with gene_length
    gene_length_fp = midas_db.get_target_layout("pangenome_tempfile", False, species_id, "vsearch", "genes.len")
    dict_of_gene_length = scan_gene_length(gene_length_fp) # <gene_id:gene_length>

    # collect SGC mapfile from all genomes of species_id
    mapfiles_by_genomes = midas_db.fetch_files("marker_genes_map", [species_id], rep_only=False)
    cat_files(mapfiles_by_genomes.values(), "mapfile", 20)
    gene_to_marker = scan_mapfile("mapfile") # <gene_id, marker_id>

    # generate genes_info.tsv
    augment_gene_info(centroid_info, gene_to_marker, dict_of_gene_length, "genes_info.tsv", )

    copy_tasks = [
        ("cdhit/genes.ffn", localpath(midas_db, species_id, "genes.ffn")),
        ("mapfile", midas_db.get_target_layout("pangenome_marker_map", False, species_id)),
        ("genes_info.tsv", localpath(midas_db, species_id, "genes_info.tsv")),
        (f"cdhit/centroids.{max_percent}.ffn", midas_db.get_target_layout("pangenome_centroids", False, species_id))
    ]

    if args.scratch_dir != ".":
        for src in flatten(cluster_files.values())[2:]:
            copy_tasks.append((src, localpath(midas_db, species_id, f"temp/{src}")))

    multithreading_map(copy_star, copy_tasks, args.num_threads)

    if not args.debug:
        command("rm -f mapfile")


def register_args(main_func):
    subparser = add_subcommand('recluster_centroids', main_func, help='Recluster the tidied centroids to operational gene family cluster levels')
    subparser.add_argument('-s',
                           '--species',
                           dest='species',
                           required=False,
                           help="species[,species...] whose pangenome(s) to build;  alternatively, species slice in format idx:modulus, e.g. 1:30, meaning recluster species whose ids are 1 mod 30; or, the special keyword 'all' meaning all species")
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
    return main_func


@register_args
def main(args):
    tsprint(f"Executing midas2 subcommand {args.subcommand}.") # with args {vars(args)}.
    recluster_centroid(args)
