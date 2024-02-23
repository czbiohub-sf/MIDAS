#!/usr/bin/env python3
import os
import sys

from midas2.common.argparser import add_subcommand
from midas2.common.utils import tsprint, retry, command, multithreading_map, find_files, upload, num_physical_cores, pythonpath, split, upload_star
from midas2.common.utilities import decode_genomes_arg, parse_gff_to_tsv
from midas2.models.midasdb import MIDAS_DB
from midas2.params.inputs import MARKER_FILE_EXTS, MIDASDB_NAMES


@retry
def find_files_with_retry(f):
    return find_files(f)


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
        local_file = midas_db.get_target_layout("annotation_genes", False, species_id, genome_id)

        msg = f"Building gene features for genome {genome_id} from species {species_id}."
        if args.upload and find_files_with_retry(dest_file):
            if not args.force:
                tsprint(f"Destination {dest_file} for genome {genome_id} gene features already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Building", "Rebuilding")
        if not args.upload and os.path.exists(local_file):
            if not args.force:
                tsprint(f"Destination {local_file} for genome {genome_id} gene features already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Building", "Rebuilding")

        tsprint(msg)
        last_dest = midas_db.get_target_layout("gene_features_log", True, species_id, genome_id)
        worker_log = midas_db.get_target_layout("gene_features_log", False, species_id, genome_id)
        worker_subdir = os.path.dirname(worker_log)

        if not args.debug:
            command(f"rm -rf {worker_subdir}")
        if not os.path.isdir(worker_subdir):
            command(f"mkdir -p {worker_subdir}")

        # Recurisve call via subcommand.  Use subdir, redirect logs.
        subcmd_str = f"--zzz_worker_mode --midasdb_name {args.midasdb_name} --midasdb_dir {os.path.abspath(args.midasdb_dir)} {'--debug' if args.debug else ''} {'--upload' if args.upload else ''}"
        worker_cmd = f"cd {worker_subdir}; PYTHONPATH={pythonpath()} {sys.executable} -m midas2 build_midasdb --generate_gene_feature --genomes {genome_id} {subcmd_str} &>> {worker_log}"
        with open(f"{worker_log}", "w") as slog:
            slog.write(msg + "\n")
            slog.write(worker_cmd + "\n")

        try:
            command(worker_cmd)
        finally:
            if args.upload:
                upload(f"{worker_log}", last_dest, check=False)
            if not args.debug:
                command(f"rm -rf {worker_subdir}", check=False)

    genome_id_list = decode_genomes_arg(args, species_for_genome)
    multithreading_map(genome_work, genome_id_list, args.num_threads)


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

    if args.upload:
        dest_file = midas_db.get_target_layout("annotation_genes", True, species_id, genome_id)
        upload(output_genes, dest_file)


def build_markerdb(args):
    """ Collate marker genes of repgenomes into phyeco.fa and phyeco.map """
    midas_db = MIDAS_DB(os.path.abspath(args.midasdb_dir), args.midasdb_name, num_cores=num_physical_cores)
    species = midas_db.uhgg.species

    last_dest = midas_db.get_target_layout("build_markerdb_log", remote=True)
    build_markerdb_log = midas_db.get_target_layout("build_markerdb_log", remote=False)

    build_marker_subdir = os.path.dirname(build_markerdb_log)
    if not args.debug:
        command(f"rm -rf {build_marker_subdir}")
    if not os.path.isdir(build_marker_subdir):
        command(f"mkdir -p {build_marker_subdir}")

    msg = "Collating marker genes sequences."
    if args.upload and find_files_with_retry(last_dest):
        if not args.force:
            tsprint(f"Destination {last_dest} already exists.  Specify --force to overwrite.")
            return
        msg = msg.replace(msg.split(" ", maxsplit=1)[0], "Re-" + msg.split(" ", maxsplit=1)[0])
    if not args.upload and os.path.exists(build_markerdb_log):
        if not args.force:
            tsprint(f"Destination {build_markerdb_log} already exists.  Specify --force to overwrite.")
            return
        msg = msg.replace(msg.split(" ", maxsplit=1)[0], "Re-" + msg.split(" ", maxsplit=1)[0])

    tsprint(msg)
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
    if args.upload:
        upload_tasks = list(zip(midas_db.get_target_layout("marker_db", False), midas_db.get_target_layout("marker_db", True)))
        multithreading_map(upload_star, upload_tasks, 4)
        upload(build_markerdb_log, last_dest, check=False) # Upload LOG file last as a checkout point


def register_args(main_func):
    subparser = add_subcommand('build_midasdb', main_func, help='Generate variety of MIDAS DB related files desired by MIDAS')
    subparser.add_argument('--genomes',
                           dest='genomes',
                           required=False,
                           help="genome[,genome...] to import;  alternatively, slice in format idx:modulus, e.g. 1:30, meaning build midasdb for genomes whose ids are 1 mod 30; or, the special keyword 'all' meaning all genomes")
    subparser.add_argument('--species',
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
    subparser.add_argument('--generate_gene_feature',
                           action='store_true',
                           default=False,
                           help="Generate gene features for each genomes")
    subparser.add_argument('--build_markerdb',
                           action='store_true',
                           default=False,
                           help="Build SGC marker db.")
    subparser.add_argument('--upload',
                           action='store_true',
                           default=False,
                           help="Upload built files to AWS S3")
    subparser.add_argument('-t',
                           '--num_threads',
                           dest='num_threads',
                           type=int,
                           default=num_physical_cores,
                           help="Number of threads")
    subparser.add_argument('--scratch_dir',
                           dest='scratch_dir',
                           type=str,
                           default=".",
                           help="Absolute path to scratch directory for fast I/O.")
    return main_func


@register_args
def main(args):
    tsprint(f"Executing midas2 subcommand {args.subcommand}.") # with args {vars(args)}.
    if args.generate_gene_feature:
        generate_gene_feature(args)
    if args.build_markerdb:
        build_markerdb(args)
