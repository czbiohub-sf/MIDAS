import os
import sys
from itertools import chain
from hashlib import md5
import Bio.SeqIO
from midas2.common.argparser import add_subcommand
from midas2.common.utils import tsprint, InputStream, retry, command, multithreading_map, find_files, upload, pythonpath, num_physical_cores
from midas2.common.utilities import decode_species_arg, decode_genomes_arg
from midas2.models.midasdb import MIDAS_DB
from midas2.params.inputs import MIDASDB_NAMES


CONCURRENT_GENOME_IMPORTS = 20


@retry
def find_files_with_retry(f):
    return find_files(f)


# 1. Occasional failures in aws s3 cp require a retry.
# 2. In future, for really large numbers of genomes, we may prefer a separate wave of retries for all first-attempt failures.
# 3. The Bio.SeqIO.parse() code is CPU-bound and thus it's best to run this function in a separate process for every genome.
@retry
def clean_genome(genome_id, raw_genome_fp):
    output_genome = f"{genome_id}.fasta"
    with open(output_genome, 'w') as o_genome, \
         InputStream(raw_genome_fp, check_path=False) as genome:
        for sn, rec in enumerate(Bio.SeqIO.parse(genome, 'fasta')):
            contig_seq = str(rec.seq).upper()
            contig_len = len(contig_seq)
            # We no longer include contigs shorter than 1000 bps into MIDASDB
            if contig_len < 1000:
                continue
            contig_hash = md5(contig_seq.encode('utf-8')).hexdigest()[-6:]
            new_contig_id = f"{genome_id}_C{sn}_L{contig_len/1000.0:3.1f}k_H{contig_hash}"
            o_genome.write(f">{new_contig_id}\n{contig_seq}\n")
    return output_genome


def import_genome(args):
    if args.zzz_worker_mode:
        import_genome_worker(args)
    else:
        import_genome_master(args)


def import_genome_master(args):

    # Fetch table of contents from s3.
    # This will be read separately by each species build subcommand, so we make a local copy.
    midas_db = MIDAS_DB(os.path.abspath(args.midasdb_dir), args.midasdb_name)
    species_for_genome = midas_db.uhgg.genomes

    def genome_work(genome_id):
        assert genome_id in species_for_genome, f"Genome {genome_id} is not in the database."
        species_id = species_for_genome[genome_id]

        dest_file = midas_db.get_target_layout("imported_genome", True, species_id, genome_id, "fasta")
        local_file = midas_db.get_target_layout("imported_genome", False, species_id, genome_id, "fasta")

        msg = f"Importing genome {genome_id} from species {species_id}."
        if args.upload and find_files_with_retry(dest_file):
            if not args.force:
                tsprint(f"Destination {dest_file} for genome {genome_id} annotations already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Importing", "Reimporting")

        if not args.upload and os.path.exists(local_file):
            if not args.force:
                tsprint(f"Destination {local_file} for genome {genome_id} annotations already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Importing", "Reimporting")

        tsprint(msg)
        last_dest = midas_db.get_target_layout("import_log", True, species_id, genome_id)
        local_dest = midas_db.get_target_layout("import_log", False, species_id, genome_id)
        local_dir = os.path.dirname(os.path.dirname(local_dest))
        command(f"mkdir -p {local_dir}")

        worker_log = os.path.basename(local_dest)
        worker_subdir = os.path.dirname(local_dest) if args.scratch_dir == "." else f"{args.scratch_dir}/import/{genome_id}"
        worker_log = f"{worker_subdir}/{worker_log}"

        if not args.debug:
            command(f"rm -rf {worker_subdir}")
        if not os.path.isdir(worker_subdir):
            command(f"mkdir -p {worker_subdir}")

        # Recurisve call via subcommand.  Use subdir, redirect logs.
        subcmd_str = f"--zzz_worker_mode --midasdb_name {args.midasdb_name} --midasdb_dir {os.path.abspath(args.midasdb_dir)} {'--debug' if args.debug else ''} {'--upload' if args.upload else ''}"
        worker_cmd = f"cd {worker_subdir}; PYTHONPATH={pythonpath()} {sys.executable} -m midas2 import_genome --genome {genome_id} {subcmd_str} &>> {worker_log}"
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
                command(f"cp -r {worker_subdir} {local_dir}")
            if not args.debug:
                command(f"rm -rf {worker_subdir}", check=False)

    if args.genomes:
        genome_id_list = decode_genomes_arg(args, species_for_genome)

    if args.species:
        species = midas_db.uhgg.species
        species_id_list = decode_species_arg(args, species)
        genome_id_list = list(chain.from_iterable([list(species[spid].keys()) for spid in species_id_list]))

    multithreading_map(genome_work, genome_id_list, num_threads=args.num_threads)


def import_genome_worker(args):

    violation = "Please do not call import_genome_worker directly.  Violation"
    assert args.zzz_worker_mode, f"{violation}:  Missing --zzz_worker_mode arg."

    midas_db = MIDAS_DB(args.midasdb_dir, args.midasdb_name)
    species_for_genome = midas_db.uhgg.genomes

    genome_id = args.genomes
    species_id = species_for_genome[genome_id]

    raw_genome_fp = midas_db.fetch_file("raw_genome", species_id, genome_id, "fa")
    cleaned = clean_genome(genome_id, raw_genome_fp)

    if args.upload:
        upload(cleaned, midas_db.get_target_layout("imported_genome", True, species_id, genome_id, "fasta"))


def register_args(main_func):
    subparser = add_subcommand('import_genome', main_func, help='import selected genomes from UHGG')
    subparser.add_argument('--genomes',
                           dest='genomes',
                           required=False,
                           help="genome[,genome...] to import;  alternatively, slice in format idx:modulus, e.g. 1:30, meaning import genomes whose ids are 1 mod 30; or, the special keyword 'all' meaning all genomes")
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
    subparser.add_argument('--scratch_dir',
                           dest='scratch_dir',
                           type=str,
                           default=".",
                           help="Absolute path to scratch directory for fast I/O.")
    return main_func


@register_args
def main(args):
    tsprint(f"Executing midas2 subcommand {args.subcommand} with args {vars(args)}.")
    import_genome(args)
