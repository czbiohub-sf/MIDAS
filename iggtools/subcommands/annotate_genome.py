#!/usr/bin/env python3
import os
import sys
from multiprocessing import Semaphore
from itertools import chain
from iggtools.common.argparser import add_subcommand, SUPPRESS
from iggtools.common.utils import tsprint, retry, command, multithreading_map, find_files, upload, pythonpath, upload_star, download_reference, drop_lz4
from iggtools.models.uhgg import UHGG, get_uhgg_layout, destpath, unified_genome_id
from iggtools.params import outputs
from iggtools.subcommands.import_uhgg import decode_genomes_arg
from iggtools.subcommands.build_pangenome import decode_species_arg


CONCURRENT_PROKKA_RUNS = Semaphore(6)


@retry
def find_files_with_retry(f):
    return find_files(f)


# 1. Occasional failures in aws s3 cp require a retry.
@retry
def download_genome(genome_id, cleaned_genome):
    command(f"rm -f {genome_id}.fasta")
    command(f"aws s3 cp --only-show-errors {cleaned_genome} - | lz4 -dc > {genome_id}.fasta")


def prokka_gene_annotation(genome_id, species_id):
    # Prokka will crash if installed <6 months ago.  It's a feature.  See tbl2asn.
    cleaned_genome = destpath(get_uhgg_layout(species_id, "fna", genome_id)["imported_genome_file"])
    ugid = unified_genome_id(genome_id)

    download_genome(genome_id, cleaned_genome)

    subdir = "prokka_dir"
    command(f"rm -rf {subdir}")

    output_files = [
        f"{genome_id}.faa",
        f"{genome_id}.ffn",
        f"{genome_id}.fna",
        f"{genome_id}.gff",
        f"{genome_id}.tsv"
    ]
    command(f"prokka --kingdom Bacteria --outdir {subdir} --cpus 8 --prefix {genome_id} --locustag {ugid} --compliant {genome_id}.fasta")
    for o in output_files:
        command(f"mv {subdir}/{o} .")

    return output_files


def annotate_genome(args):
    if args.zzz_worker_toc:
        annotate_genome_worker(args)
    else:
        annotate_genome_master(args)


def annotate_genome_master(args):

    # Fetch table of contents from s3.
    # This will be read separately by each species build subcommand, so we make a local copy.
    local_toc = download_reference(outputs.genomes)

    db = UHGG(local_toc)
    species_for_genome = db.genomes

    def genome_work(genome_id):
        assert genome_id in species_for_genome, f"Genome {genome_id} is not in the database."
        species_id = species_for_genome[genome_id]

        dest_file = destpath(get_uhgg_layout(species_id, "fna", genome_id)["annotation_file"])
        msg = f"Annotating genome {genome_id} from species {species_id}."
        if find_files_with_retry(dest_file):
            if not args.force:
                tsprint(f"Destination {dest_file} for genome {genome_id} annotations already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Importing", "Reimporting")


        with CONCURRENT_PROKKA_RUNS:

            tsprint(msg)
            logfile = get_uhgg_layout(species_id, "", genome_id)["annotation_log"]
            worker_log = os.path.basename(logfile)
            worker_subdir = f"{species_id}__{genome_id}"
            if not args.debug:
                command(f"rm -rf {worker_subdir}")
            if not os.path.isdir(worker_subdir):
                command(f"mkdir {worker_subdir}")
            # Recurisve call via subcommand.  Use subdir, redirect logs.
            worker_cmd = f"cd {worker_subdir}; PYTHONPATH={pythonpath()} {sys.executable} -m iggtools annotate_genome --genome {genome_id} --zzz_worker_mode --zzz_worker_toc {os.path.abspath(local_toc)} {'--debug' if args.debug else ''} &>> {worker_log}"
            with open(f"{worker_subdir}/{worker_log}", "w") as slog:
                slog.write(msg + "\n")
                slog.write(worker_cmd + "\n")
            try:
                command(worker_cmd)
            finally:
                # Cleanup should not raise exceptions of its own, so as not to interfere with any
                # prior exceptions that may be more informative.  Hence check=False.
                upload(f"{worker_subdir}/{worker_log}", destpath(logfile), check=False)
                if not args.debug:
                    command(f"rm -rf {worker_subdir}", check=False)

    if args.genomes:
        genome_id_list = decode_genomes_arg(args, species_for_genome)
    if args.species:
        species = db.species
        species_id_list = decode_species_arg(args, species)
        genome_id_list = list(chain.from_iterable([list(species[spid].keys()) for spid in species_id_list]))

    multithreading_map(genome_work, genome_id_list, num_threads=20)


def annotate_genome_worker(args):
    """
    https://github.com/czbiohub/iggtools/wiki
    """

    violation = "Please do not call annotate_genome_worker directly.  Violation"
    assert args.zzz_worker_mode, f"{violation}:  Missing --zzz_worker_mode arg."
    assert os.path.isfile(args.zzz_worker_toc), f"{violation}: File does not exist: {args.zzz_worker_toc}"

    db = UHGG(args.zzz_worker_toc)
    species_for_genome = db.genomes

    genome_id = args.genomes
    species_id = species_for_genome[genome_id]

    dest_file = destpath(get_uhgg_layout(species_id, "fna", genome_id)["annotation_file"])
    last_output = os.path.basename(dest_file)

    output_files = prokka_gene_annotation(genome_id, species_id)

    upload_tasks = []
    for o in output_files:
        olz = o + ".lz4"
        otype = o.rsplit(".")[-1]
        if olz != last_output:
            upload_tasks.append((o, destpath(get_uhgg_layout(species_id, otype, genome_id)["annotation_file"])))

    command(f"aws s3 rm --recursive {dest_file.rsplit('/', 1)[0]}")
    multithreading_map(upload_star, upload_tasks)

    # Upload this last because it indicates all other work has succeeded.
    upload(drop_lz4(last_output), dest_file)


def register_args(main_func):
    subparser = add_subcommand('annotate_genome', main_func, help='annotate selected genomes with prokka')
    subparser.add_argument('--genomes',
                           dest='genomes',
                           required=False,
                           help="genome[,genome...] to import;  alternatively, slice in format idx:modulus, e.g. 1:30, meaning annotate genomes whose ids are 1 mod 30; or, the special keyword 'all' meaning all genomes")
    subparser.add_argument('--species',
                           dest='species',
                           required=False,
                           help="species[,species...] whose pangenome(s) to build;  alternatively, species slice in format idx:modulus, e.g. 1:30, meaning build species whose ids are 1 mod 30; or, the special keyword 'all' meaning all species")
    subparser.add_argument('--zzz_worker_toc',
                           dest='zzz_worker_toc',
                           required=False,
                           help=SUPPRESS) # "reserved to pass table of contents from master to worker"
    return main_func


@register_args
def main(args):
    tsprint(f"Executing iggtools subcommand {args.subcommand} with args {vars(args)}.")
    annotate_genome(args)
