#!/usr/bin/env python3
import os
import sys
from multiprocessing import Semaphore
import gffutils

from midas2.common.argparser import add_subcommand
from midas2.common.utils import tsprint, OutputStream, retry, command, multithreading_map, find_files, upload, num_physical_cores, pythonpath, cat_files, split, upload_star
from midas2.common.utilities import scan_mapfile, scan_gene_info, scan_gene_length
from midas2.common.utilities import decode_species_arg, decode_genomes_arg
from midas2.models.midasdb import MIDAS_DB
from midas2.params.schemas import CLUSTER_INFO_SCHEMA
from midas2.params.inputs import MARKER_FILE_EXTS, MIDASDB_NAMES


# Up to this many concurrent species builds.
CONCURRENT_BUILDS = Semaphore(num_physical_cores)


@retry
def find_files_with_retry(f):
    return find_files(f)


def parse_gff_to_tsv(gff3_file, genes_file):
    """ Convert GFF3 features format into genes.feature """
    db = gffutils.create_db(gff3_file, f"{gff3_file}.db")
    with OutputStream(genes_file) as stream:
        stream.write("\t".join(["gene_id", "contig_id", "start", "end", "strand", "gene_type"]) + "\n")
        for feature in db.all_features():
            if feature.source == "prokka":
                continue
            if "ID" not in feature.attributes: #CRISPR
                continue
            seqid = feature.seqid
            start = feature.start
            stop = feature.stop
            strand = feature.strand
            gene_id = feature.attributes['ID'][0]
            locus_tag = feature.attributes['locus_tag'][0]
            assert gene_id == locus_tag
            gene_type = feature.featuretype
            stream.write("\t".join([gene_id, seqid, str(start), str(stop), strand, gene_type]) + "\n")
    return True


def generate_gene_feature(args):
    if args.zzz_worker_mode:
        generate_gene_feature_worker(args)
    else:
        generate_gene_feature_master(args)


def generate_gene_feature_master(args):

    midas_db = MIDAS_DB(os.path.abspath(args.midasdb_dir), args.midasdb_name)
    species_for_genome = midas_db.uhgg.genomes

    def genome_work(genome_id):
        assert genome_id in species_for_genome, f"Genome {genome_id} is not in the database."
        species_id = species_for_genome[genome_id]

        dest_file = midas_db.get_target_layout("annotation_genes", True, species_id, genome_id)
        msg = f"Builing gene features for genome {genome_id} from species {species_id}."
        if find_files_with_retry(dest_file):
            if not args.force:
                tsprint(f"Destination {dest_file} for genome {genome_id} gene features already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Importing", "Reimporting")

        with CONCURRENT_BUILDS:
            tsprint(msg)
            last_dest = midas_db.get_target_layout("gene_features_log", True, species_id, genome_id)

            worker_log = midas_db.get_target_layout("gene_features_log", False, species_id, genome_id)
            worker_subdir = os.path.dirname(worker_log)

            if not args.debug:
                command(f"rm -rf {worker_subdir}")
            if not os.path.isdir(worker_subdir):
                command(f"mkdir -p {worker_subdir}")

            # Recurisve call via subcommand.  Use subdir, redirect logs.
            worker_cmd = f"cd {worker_subdir}; PYTHONPATH={pythonpath()} {sys.executable} -m midas2 build_midasdb --generate_gene_feature --genomes {genome_id} --midasdb_name {args.midasdb_name} --midasdb_dir {os.path.abspath(args.midasdb_dir)} --zzz_worker_mode {'--debug' if args.debug else ''} &>> {worker_log}"
            with open(f"{worker_log}", "w") as slog:
                slog.write(msg + "\n")
                slog.write(worker_cmd + "\n")

            try:
                command(worker_cmd)
            finally:
                if not args.debug:
                    upload(f"{worker_log}", last_dest, check=False)
                    command(f"rm -rf {worker_subdir}", check=False)

    genome_id_list = decode_genomes_arg(args, species_for_genome)
    multithreading_map(genome_work, genome_id_list, num_physical_cores)


def generate_gene_feature_worker(args):
    """
    https://github.com/czbiohub/MIDAS2.0/wiki/MIDAS-DB
    """

    violation = "Please do not call generate_gene_feature_worker directly.  Violation"
    assert args.zzz_worker_mode, f"{violation}:  Missing --zzz_worker_mode arg."

    genome_id = args.genomes
    midas_db = MIDAS_DB(args.midasdb_dir, args.midasdb_name)

    species_for_genome = midas_db.uhgg.genomes
    species_id = species_for_genome[genome_id]

    input_annotation = midas_db.fetch_file("annotation_gff", species_id, genome_id)
    output_genes = midas_db.get_target_layout("annotation_genes", False, species_id, genome_id)
    assert parse_gff_to_tsv(input_annotation, output_genes)

    if not args.debug:
        dest_file = midas_db.get_target_layout("annotation_genes", True, species_id, genome_id)
        upload(output_genes, dest_file)


def generate_cluster_info(args):
    if args.zzz_worker_mode: #"zzz_worker_mode" in args and args.zzz_worker_mode:
        generate_cluster_info_worker(args)
    else:
        generate_cluster_info_master(args)


def generate_cluster_info_master(args):

    midas_db = MIDAS_DB(os.path.abspath(args.midasdb_dir), args.midasdb_name)
    species = midas_db.uhgg.species

    def species_work(species_id):
        assert species_id in species, f"Species {species_id} is not in the database."
        species_genomes = species[species_id]

        dest_file = midas_db.get_target_layout("pangenome_cluster_info", True, species_id)
        msg = f"Computing cluster_info for species {species_id} with {len(species_genomes)} total genomes."
        if find_files_with_retry(dest_file):
            if not args.force:
                tsprint(f"Destination {dest_file} for species {species_id} cluster_info already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Computing", "Recomputing")

        with CONCURRENT_BUILDS:
            tsprint(msg)
            last_dest = midas_db.get_target_layout("cluster_info_log", True, species_id)

            worker_log = midas_db.get_target_layout("cluster_info_log", False, species_id)
            worker_subdir = os.path.dirname(worker_log)

            if not args.debug:
                command(f"rm -rf {worker_subdir}")
            if not os.path.isdir(worker_subdir):
                command(f"mkdir -p {worker_subdir}")

            # Recurisve call via subcommand.  Use subdir, redirect logs.
            worker_cmd = f"cd {worker_subdir}; PYTHONPATH={pythonpath()} {sys.executable} -m midas2 build_midasdb --generate_cluster_info --species {species_id} --midasdb_name {args.midasdb_name} --midasdb_dir {os.path.abspath(args.midasdb_dir)} --zzz_worker_mode {'--debug' if args.debug else ''} &>> {worker_log}"
            with open(f"{worker_log}", "w") as slog:
                slog.write(msg + "\n")
                slog.write(worker_cmd + "\n")

            try:
                command(worker_cmd)
            finally:
                if not args.debug:
                    upload(f"{worker_log}", last_dest, check=False)
                    command(f"rm -rf {worker_subdir}", check=False)

    species_id_list = decode_species_arg(args, species)
    multithreading_map(species_work, species_id_list, num_physical_cores)


def generate_cluster_info_worker(args):

    violation = "Please do not call generate_cluster_info_worker directly.  Violation"
    assert args.zzz_worker_mode, f"{violation}:  Missing --zzz_worker_mode arg."

    species_id = args.species
    midas_db = MIDAS_DB(args.midasdb_dir, args.midasdb_name)
    species = midas_db.uhgg.species
    assert species_id in species, f"{violation}: Species {species_id} is not in the database."

    # Collect mapfile from all genomes of species_id
    list_of_genomes = species[species_id].keys()
    list_of_mapfiles = [midas_db.fetch_file("marker_genes_map", species_id, genome_id) for genome_id in list_of_genomes]
    species_mapfile = f"mapfile"
    cat_files(list_of_mapfiles, species_mapfile, 20)

    # <gene_id, marker_id>
    dict_of_markers = scan_mapfile(species_mapfile)

    # <centroid_99> as the key
    gene_info_fp = midas_db.fetch_file("pangenome_genes_info", species_id)
    dict_of_centroids = scan_gene_info(gene_info_fp)

    # <gene_id:gene_length>
    gene_length_file = midas_db.fetch_file("pangenome_genes_len", species_id)
    dict_of_gene_length = scan_gene_length(gene_length_file)

    dest_file = midas_db.get_target_layout("pangenome_cluster_info", True, species_id)
    local_file = midas_db.get_target_layout("pangenome_cluster_info", False, species_id)

    with OutputStream(local_file) as stream:
        stream.write("\t".join(CLUSTER_INFO_SCHEMA.keys()) + "\n")
        for r in dict_of_centroids.values():
            centroid_99 = r["centroid_99"]

            marker_id = dict_of_markers[centroid_99] if centroid_99 in dict_of_markers else ""
            gene_len = dict_of_gene_length[centroid_99]

            val = list(r.values()) + [gene_len, marker_id]
            stream.write("\t".join(map(str, val)) + "\n")

    if not args.debug:
        command(f"aws s3 rm {dest_file}")
        upload(local_file, dest_file, check=False)


def build_markerdb(args):
    """ Collate marker genes of repgenomes into phyeco.fa and phyeco.map """
    midas_db = MIDAS_DB(os.path.abspath(args.midasdb_dir), args.midasdb_name, num_cores=num_physical_cores)
    species = midas_db.uhgg.species

    build_markerdb_log_s3 = midas_db.get_target_layout("build_markerdb_log", remote=True)
    msg = f"Collating marker genes sequences."
    if find_files_with_retry(build_markerdb_log_s3):
        if not args.force:
            tsprint(f"Destination {build_markerdb_log_s3} already exists.  Specify --force to overwrite.")
            return
        msg = msg.replace(msg.split(" ")[0], "Re-" + msg.split(" ")[0])
    tsprint(msg)

    build_markerdb_log = midas_db.get_target_layout("build_markerdb_log", remote=False)
    build_marker_subdir = os.path.dirname(build_markerdb_log)

    if not args.debug:
        command(f"rm -rf {build_marker_subdir}")
    if not os.path.isdir(build_marker_subdir):
        command(f"mkdir -p {build_marker_subdir}")

    with open(build_markerdb_log, "w") as slog:
        slog.write(msg + "\n")

    # Fetch marker genes fasta-file and map-file from s3
    list_of_repgenome_marker_seqs = midas_db.fetch_files("marker_genes_seq", species.keys())
    list_of_repgenome_marker_maps = midas_db.fetch_files("marker_genes_map", species.keys())

    # Collate to phyeco.fa and phyeco.map
    phyeco_files = dict(zip(MARKER_FILE_EXTS, midas_db.get_target_layout("marker_db", False)))
    phyeco_seqs = phyeco_files["fa"]
    for _ in split(list_of_repgenome_marker_seqs.values(), 20):
        command("cat " + " ".join(_) + f" >> {phyeco_seqs}")

    phyeco_maps = phyeco_files["map"]
    for _ in split(list_of_repgenome_marker_maps.values(), 20):
        command("cat " + " ".join(_) + f" >> {phyeco_maps}")

    # Build hs-blastn index for the collated phyeco sequences
    cmd_index = f"hs-blastn index {phyeco_seqs} &>> {build_markerdb_log}"
    with open(f"{build_markerdb_log}", "a") as slog:
        slog.write(cmd_index + "\n")
    command(cmd_index)

    # Upload generated fasta and index files
    if not args.debug:
        upload_tasks = list(zip(midas_db.get_target_layout("marker_db", False), midas_db.get_target_layout("marker_db", True)))
        multithreading_map(upload_star, upload_tasks)
        upload(build_markerdb_log, build_markerdb_log_s3, check=False) # Upload LOG file last as a checkout point


def register_args(main_func):
    subparser = add_subcommand('build_midasdb', main_func, help='Generate variety of MIDAS DB related files desired by MIDAS')
    subparser.add_argument('--genomes',
                           dest='genomes',
                           required=False,
                           help="genome[,genome...] to import;  alternatively, slice in format idx:modulus, e.g. 1:30, meaning annotate genomes whose ids are 1 mod 30; or, the special keyword 'all' meaning all genomes")
    subparser.add_argument('--species',
                           dest='species',
                           required=False,
                           help="species[,species...] whose pangenome(s) to build;  alternatively, species slice in format idx:modulus, e.g. 1:30, meaning build species whose ids are 1 mod 30; or, the special keyword 'all' meaning all species")
    subparser.add_argument('--midasdb_name',
                           dest='midasdb_name',
                           type=str,
                           default="uhgg",
                           choices=MIDASDB_NAMES,
                           help=f"MIDAS Database name.")
    subparser.add_argument('--midasdb_dir',
                           dest='midasdb_dir',
                           type=str,
                           default=".",
                           help=f"MIDAS Database name.")
    subparser.add_argument('--generate_gene_feature',
                           action='store_true',
                           default=False,
                           help=f"Generate gene features for each genomes")
    subparser.add_argument('--generate_cluster_info',
                           action='store_true',
                           default=False,
                           help=f"Generate cluster_info.txt used in chunk design and merge_genes.")
    subparser.add_argument('--build_markerdb',
                           action='store_true',
                           default=False,
                           help=f"Generate cluster_info.txt used in chunk design and merge_genes.")
    return main_func


@register_args
def main(args):
    tsprint(f"Executing midas2 subcommand {args.subcommand} with args {vars(args)}.")
    if args.generate_gene_feature:
        generate_gene_feature(args)
    if args.generate_cluster_info:
        generate_cluster_info(args)
    if args.build_markerdb:
        build_markerdb(args)
