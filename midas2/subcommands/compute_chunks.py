#!/usr/bin/env python3
import os
import sys
import json

from midas2.common.argparser import add_subcommand
from midas2.common.utils import tsprint, num_physical_cores, retry, find_files, pythonpath, upload, OutputStream, command, multithreading_map
from midas2.common.utilities import decode_species_arg
from midas2.models.midasdb import MIDAS_DB
from midas2.models.species import design_run_snps_chunks, design_merge_snps_chunks
from midas2.params.inputs import MIDASDB_NAMES


DEFAULT_CHUNK_SIZE = 100000


@retry
def find_files_with_retry(f):
    return find_files(f)


def get_dest_filename(chunk_type, species_id, genome_id=""):
    if chunk_type == "run_snps":
        dest_filename = "chunks_sites_run"
        msg = f"Designing chunks of sites for representative genome {genome_id} from species {species_id} for RUN snps."

    if chunk_type == "merge_snps":
        dest_filename = "chunks_sites_merge"
        msg = f"Designing chunks of sites for representative genome {genome_id} from species {species_id} for MERGE snps."

    return (dest_filename, msg)


def compute_chunks(args):
    if args.zzz_worker_mode:
        compute_chunks_worker(args)
    else:
        compute_chunks_master(args)


def compute_chunks_master(args):

    midas_db = MIDAS_DB(os.path.abspath(args.midasdb_dir), args.midasdb_name)
    repgenome_for_species = midas_db.uhgg.representatives

    def species_work(species_id):
        assert species_id in repgenome_for_species, f"Species {species_id} is not in the database."
        genome_id = repgenome_for_species[species_id]

        dest_filename, msg = get_dest_filename(args.chunk_type, species_id, genome_id)
        dest_file = midas_db.get_target_layout(dest_filename, True, species_id, genome_id, args.chunk_size)
        local_file = midas_db.get_target_layout(dest_filename, False, species_id, genome_id, args.chunk_size)

        if args.upload and find_files_with_retry(dest_file):
            if not args.force:
                tsprint(f"Destination {dest_file} for species {species_id} already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Designing", "Redesigning")
        if not args.upload and os.path.exists(local_file):
            if not args.force:
                tsprint(f"Destination {local_file} for species {species_id} already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Designing", "Redesigning")

        tsprint(msg)
        local_file = midas_db.get_target_layout(dest_filename, False, species_id, genome_id, args.chunk_size)

        worker_log = f"{species_id}_{dest_filename}.log"
        worker_subdir = os.path.dirname(local_file)

        if not args.debug:
            command(f"rm -f {local_file}")
        if not os.path.isdir(worker_subdir):
            command(f"mkdir -p {worker_subdir}")

        # Recurisve call via subcommand.  Use subdir, redirect logs.
        subcmd_str = f"--zzz_worker_mode --midasdb_name {args.midasdb_name} --midasdb_dir {os.path.abspath(args.midasdb_dir)} {'--debug' if args.debug else ''} {'--upload' if args.upload else ''}"
        worker_cmd = f"cd {worker_subdir}; PYTHONPATH={pythonpath()} {sys.executable} -m midas2 compute_chunks --species {species_id} --chunk_size {args.chunk_size} --chunk_type {args.chunk_type} {subcmd_str} &>> {worker_subdir}/{worker_log}"
        with open(f"{worker_subdir}/{worker_log}", "w") as slog:
            slog.write(msg + "\n")
            slog.write(worker_cmd + "\n")
        try:
            command(worker_cmd)
        finally:
            if not args.debug:
                command(f"rm -rf {worker_subdir}/{worker_log} {local_file}", check=False)

    species_id_list = decode_species_arg(args, repgenome_for_species)
    multithreading_map(species_work, species_id_list, args.num_threads)


def compute_chunks_worker(args):

    violation = "Please do not call compute_chunks_worker directly.  Violation"
    assert args.zzz_worker_mode, f"{violation}:  Missing --zzz_worker_mode arg."

    species_id = args.species
    midas_db = MIDAS_DB(args.midasdb_dir, args.midasdb_name)

    repgenome_for_species = midas_db.uhgg.representatives
    assert species_id in repgenome_for_species, f"{violation}: Species {species_id} is not in the database."
    genome_id = repgenome_for_species[species_id]

    if args.chunk_type == "run_snps":
        contigs_fp = midas_db.fetch_file("representative_genome", species_id, genome_id)
        chunks_to_cache = design_run_snps_chunks(species_id, contigs_fp, args.chunk_size)
    if args.chunk_type == "merge_snps":
        contigs_fp = midas_db.fetch_file("representative_genome", species_id, genome_id)
        chunks_to_cache = design_merge_snps_chunks(species_id, contigs_fp, args.chunk_size)

    dest_filename, _ = get_dest_filename(args.chunk_type, species_id, genome_id)

    local_file = midas_db.get_target_layout(dest_filename, False, species_id, genome_id, args.chunk_size)
    with OutputStream(local_file) as stream:
        json.dump(chunks_to_cache, stream)

    dest_file = midas_db.get_target_layout(dest_filename, True, species_id, genome_id, args.chunk_size)
    if args.upload:
        command(f"aws s3 rm {dest_file}")
        upload(local_file, dest_file)


def register_args(main_func):
    subparser = add_subcommand('compute_chunks', main_func, help='Design chunks for SNPs or Genes for given species')
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
    subparser.add_argument('--chunk_type',
                           dest='chunk_type',
                           type=str,
                           default="run_snps",
                           choices=['run_snps', 'merge_snps', "genes"],
                           help="chunk_type of chunks to compute.")
    subparser.add_argument('--chunk_size',
                           dest='chunk_size',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_CHUNK_SIZE,
                           help=f"Number of genomic sites for the temporary chunk file  ({DEFAULT_CHUNK_SIZE})")
    subparser.add_argument('-t',
                           '--num_threads',
                           dest='num_threads',
                           type=int,
                           default=num_physical_cores,
                           help="Number of threads")
    subparser.add_argument('--upload',
                           action='store_true',
                           default=False,
                           help="Upload built files to AWS S3")
    return main_func


@register_args
def main(args):
    tsprint(f"Executing midas2 subcommand {args.subcommand}.") # with args {vars(args)}.
    compute_chunks(args)
