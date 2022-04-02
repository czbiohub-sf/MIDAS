#!/usr/bin/env python3
import os
import sys
from multiprocessing import Semaphore
from itertools import chain
from midas2.common.argparser import add_subcommand
from midas2.common.utils import tsprint, retry, command, multithreading_map, drop_lz4, find_files, upload, pythonpath, upload_star, num_physical_cores
from midas2.common.utilities import decode_species_arg, decode_genomes_arg
from midas2.models.uhgg import unified_genome_id
from midas2.models.midasdb import MIDAS_DB
from midas2.params.inputs import MIDASDB_NAMES


CONCURRENT_PROKKA_RUNS = Semaphore(6)


@retry
def find_files_with_retry(f):
    return find_files(f)


# 1. Occasional failures in aws s3 cp require a retry.
@retry
def rename_genome(genome_id, cleaned_genome):
    command(f"rm -f {genome_id}.fasta")
    command(f"aws s3 cp --only-show-errors {cleaned_genome} - | lz4 -dc > {genome_id}.fasta")


def run_prokka(genome_id, cleaned_genome):
    # Prokka will crash if installed <6 months ago.  It's a feature.  See tbl2asn.

    ugid = unified_genome_id(genome_id) #<---- this should be obsolete for the next generation of database
    rename_genome(genome_id, cleaned_genome)

    subdir = "prokka_dir"
    command(f"rm -rf {subdir}")

    output_files = [
        f"{genome_id}.faa",
        f"{genome_id}.ffn",
        f"{genome_id}.fna",
        f"{genome_id}.gff",
        f"{genome_id}.tsv"
    ]
    command(f"prokka --kingdom Bacteria --outdir {subdir} --cpus 8 --prefix {genome_id} --locustag {genome_id} --compliant {genome_id}.fasta") #ugid
    for o in output_files:
        command(f"mv {subdir}/{o} .")

    return output_files


def annotate_genome(args):
    if args.zzz_worker_mode:
        annotate_genome_worker(args)
    else:
        annotate_genome_master(args)


def annotate_genome_master(args):

    # Fetch table of contents from s3.
    # This will be read separately by each species build subcommand, so we make a local copy.

    midas_db = MIDAS_DB(os.path.abspath(args.midasdb_dir), args.midasdb_name)
    species_for_genome = midas_db.uhgg.genomes

    def genome_work(genome_id):
        assert genome_id in species_for_genome, f"Genome {genome_id} is not in the database."
        species_id = species_for_genome[genome_id]

        dest_file = midas_db.get_target_layout("annotation_file", True, species_id, genome_id, "fna")

        msg = f"Annotating genome {genome_id} from species {species_id}."
        if find_files_with_retry(dest_file):
            if not args.force:
                tsprint(f"Destination {dest_file} for genome {genome_id} annotations already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Importing", "Reimporting")


        with CONCURRENT_PROKKA_RUNS:
            tsprint(msg)
            last_dest = midas_db.get_target_layout("annotation_log", True, species_id, genome_id)

            worker_log = midas_db.get_target_layout("annotation_log", False, species_id, genome_id)
            worker_subdir = os.path.dirname(worker_log)

            if not args.debug:
                command(f"rm -rf {worker_subdir}")
            if not os.path.isdir(worker_subdir):
                command(f"mkdir -p {worker_subdir}")

            # Recurisve call via subcommand.  Use subdir, redirect logs.
            worker_cmd = f"cd {worker_subdir}; PYTHONPATH={pythonpath()} {sys.executable} -m midas2 annotate_genome --genome {genome_id} --zzz_worker_mode --midasdb_name {args.midasdb_name} --midasdb_dir {os.path.abspath(args.midasdb_dir)} {'--debug' if args.debug else ''} &>> {worker_log}"
            with open(f"{worker_log}", "w") as slog:
                slog.write(msg + "\n")
                slog.write(worker_cmd + "\n")

            try:
                command(worker_cmd)
            finally:
                # Cleanup should not raise exceptions of its own, so as not to interfere with any
                # prior exceptions that may be more informative.  Hence check=False.
                if not args.debug:
                    upload(f"{worker_log}", last_dest, check=False)
                    command(f"rm -rf {worker_subdir}", check=False)

    if args.genomes:
        genome_id_list = decode_genomes_arg(args, species_for_genome)

    if args.species:
        species = midas_db.uhgg.species
        species_id_list = decode_species_arg(args, species)
        genome_id_list = list(chain.from_iterable([list(species[spid].keys()) for spid in species_id_list]))

    multithreading_map(genome_work, genome_id_list, num_threads=num_physical_cores)


def annotate_genome_worker(args):
    """
    https://github.com/czbiohub/midas2/wiki
    """

    violation = "Please do not call annotate_genome_worker directly.  Violation"
    assert args.zzz_worker_mode, f"{violation}:  Missing --zzz_worker_mode arg."

    midas_db = MIDAS_DB(args.midasdb_dir, args.midasdb_name)
    species_for_genome = midas_db.uhgg.genomes

    genome_id = args.genomes
    species_id = species_for_genome[genome_id]

    dest_file = midas_db.get_target_layout("annotation_file", True, species_id, genome_id, "fna")
    last_output = drop_lz4(os.path.basename(dest_file))

    cleaned_genome_fp = midas_db.get_target_layout("imported_genome", True, species_id, genome_id, "fna")
    output_files = run_prokka(genome_id, cleaned_genome_fp)

    if not args.debug:
        upload_tasks = []
        for o in output_files:
            otype = o.rsplit(".")[-1]
            if o != last_output:
                upload_tasks.append((o, midas_db.get_target_layout("annotation_file", True, species_id, genome_id, otype)))

        command(f"aws s3 rm --recursive {dest_file.rsplit('/', 1)[0]}")
        multithreading_map(upload_star, upload_tasks)

        # Upload this last because it indicates all other work has succeeded.
        upload(last_output, dest_file)


def register_args(main_func):
    subparser = add_subcommand('annotate_genome', main_func, help='Genome annotation for specified genomes using Prokka')
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
                           help=f"Local MIDAS Database path mirroing S3.")
    return main_func


@register_args
def main(args):
    tsprint(f"Executing midas2 subcommand {args.subcommand} with args {vars(args)}.")
    annotate_genome(args)
