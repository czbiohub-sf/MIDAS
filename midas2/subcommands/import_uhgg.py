#!/usr/bin/env python3
import os
import sys
from hashlib import md5
import Bio.SeqIO
from midas2.common.argparser import add_subcommand, SUPPRESS
from midas2.common.utils import tsprint, InputStream, retry, command, multithreading_map, find_files, upload, pythonpath
from midas2.common.utilities import decode_genomes_arg
from midas2.models.uhgg import UHGG, get_uhgg_layout, unified_genome_id
from midas2.params import outputs
from midas2.params.inputs import MIDASDB_DICT


CONCURRENT_GENOME_IMPORTS = 20


@retry
def find_files_with_retry(f):
    return find_files(f)


def destpath(local_path):
    igg = MIDASDB_DICT["uhgg"]
    return os.path.join(igg, f"{local_path}.lz4")


# 1. Occasional failures in aws s3 cp require a retry.
# 2. In future, for really large numbers of genomes, we may prefer a separate wave of retries for all first-attempt failures.
# 3. The Bio.SeqIO.parse() code is CPU-bound and thus it's best to run this function in a separate process for every genome.
@retry
def clean_genome(genome_id, representative_id):
    raw_genome = get_uhgg_layout(representative_id, "fna.lz4", genome_id)["raw_genome_file"]
    output_genome = f"{genome_id}.fna"

    with open(output_genome, 'w') as o_genome, \
         InputStream(raw_genome, check_path=False) as genome:
        for sn, rec in enumerate(Bio.SeqIO.parse(genome, 'fasta')):
            contig_seq = str(rec.seq).upper()
            contig_len = len(contig_seq)
            ugid = unified_genome_id(genome_id)
            contig_hash = md5(contig_seq.encode('utf-8')).hexdigest()[-6:]
            new_contig_id = f"{ugid}_C{sn}_L{contig_len/1000.0:3.1f}k_H{contig_hash}"
            o_genome.write(f">{new_contig_id}\n{contig_seq}\n")

    return output_genome


def import_uhgg(args):
    if args.zzz_worker_toc:
        import_uhgg_worker(args)
    else:
        import_uhgg_master(args)


def import_uhgg_master(args):

    # Fetch table of contents from s3.
    # This will be read separately by each species build subcommand, so we make a local copy.
    output_genomes = outputs.genomes()
    local_toc = os.path.basename(output_genomes)
    command(f"rm -f {local_toc}")
    command(f"aws s3 cp --only-show-errors {output_genomes} {local_toc}")

    db = UHGG(local_toc)
    species_for_genome = db.genomes

    def genome_work(genome_id):
        assert genome_id in species_for_genome, f"Genome {genome_id} is not in the database."
        species_id = species_for_genome[genome_id]

        dest_file = destpath(get_uhgg_layout(species_id, "fna", genome_id)["imported_genome"])

        msg = f"Importing genome {genome_id} from species {species_id}."
        if find_files_with_retry(dest_file):
            if not args.force:
                tsprint(f"Destination {dest_file} for genome {genome_id} already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Importing", "Reimporting")

        tsprint(msg)
        logfile = get_uhgg_layout(species_id, "", genome_id)["imported_genome_log"]
        worker_log = os.path.basename(logfile)
        worker_subdir = f"{species_id}__{genome_id}"
        if not args.debug:
            command(f"rm -rf {worker_subdir}")
        if not os.path.isdir(worker_subdir):
            command(f"mkdir {worker_subdir}")
        # Recurisve call via subcommand.  Use subdir, redirect logs.
        worker_cmd = f"cd {worker_subdir}; PYTHONPATH={pythonpath()} {sys.executable} -m midas2 import_uhgg --genome {genome_id} --zzz_worker_mode --zzz_worker_toc {os.path.abspath(local_toc)} {'--debug' if args.debug else ''} &>> {worker_log}"
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

    genome_id_list = decode_genomes_arg(args, species_for_genome)
    multithreading_map(genome_work, genome_id_list, num_threads=CONCURRENT_GENOME_IMPORTS)


def import_uhgg_worker(args):
    """
    https://github.com/czbiohub/MIDAS2.0/wiki/MIDAS-DB
    """

    violation = "Please do not call import_uhgg_worker directly.  Violation"
    assert args.zzz_worker_mode, f"{violation}:  Missing --zzz_worker_mode arg."
    assert os.path.isfile(args.zzz_worker_toc), f"{violation}: File does not exist: {args.zzz_worker_toc}"

    db = UHGG(args.zzz_worker_toc)
    representatives = db.representatives
    species_for_genome = db.genomes

    genome_id = args.genomes
    species_id = species_for_genome[genome_id]
    representative_id = representatives[species_id]

    dest = destpath(get_uhgg_layout(species_id, "fna", genome_id)["imported_genome"])
    command(f"aws s3 rm --recursive {os.path.dirname(dest)}")
    cleaned = clean_genome(genome_id, representative_id)
    upload(cleaned, dest)


def register_args(main_func):
    subparser = add_subcommand('import_uhgg', main_func, help='Import selected genomes from UHGG')
    subparser.add_argument('--genomes',
                           dest='genomes',
                           required=False,
                           help="genome[,genome...] to import;  alternatively, slice in format idx:modulus, e.g. 1:30, meaning import genomes whose ids are 1 mod 30; or, the special keyword 'all' meaning all genomes")
    subparser.add_argument('--zzz_worker_toc',
                           dest='zzz_worker_toc',
                           required=False,
                           help=SUPPRESS) # "reserved to pass table of contents from master to worker"
    return main_func


@register_args
def main(args):
    tsprint(f"Executing midas2 subcommand {args.subcommand} with args {vars(args)}.")
    import_uhgg(args)
