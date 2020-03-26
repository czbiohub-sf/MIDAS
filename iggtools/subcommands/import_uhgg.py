import os
import sys
from hashlib import md5
import Bio.SeqIO
from iggtools.common.argparser import add_subcommand, SUPPRESS
from iggtools.common.utils import tsprint, InputStream, retry, command, multithreading_map, find_files, upload, pythonpath
from iggtools.models.uhgg import UHGG, imported_genome_file, raw_genome_file
from iggtools.params import inputs, outputs


CONCURRENT_GENOME_IMPORTS = 20


# Move genome id parsing and name transformations in some central place that all commands can import
def unified_genome_id(genome_id):
    return "UHGG" + genome_id.replace("GUT_GENOME", "")


# 1. Occasional failures in aws s3 cp require a retry.
# 2. In future, for really large numbers of genomes, we may prefer a separate wave of retries for all first-attempt failures.
# 3. The Bio.SeqIO.parse() code is CPU-bound and thus it's best to run this function in a separate process for every genome.
@retry
def clean_genome(genome_id, representative_id):
    raw_genome = raw_genome_file(genome_id, representative_id)
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
    if args.zzz_slave_toc:
        import_uhgg_slave(args)
    else:
        import_uhgg_master(args)


@retry
def find_files_with_retry(f):
    return find_files(f)


def decode_genomes_arg(args, genomes):
    selected_genomes = set()
    try:  # pylint: disable=too-many-nested-blocks
        if args.genomes.upper() == "ALL":
            selected_genomes = set(genomes)
        else:
            for g in args.genomes.split(","):
                if ":" not in g:
                    selected_genomes.add(g)
                else:
                    i, n = g.split(":")
                    i = int(i)
                    n = int(n)
                    assert 0 <= i < n, f"Genome class and modulus make no sense: {i}, {n}"
                    for gid in genomes:
                        gid_int = int(gid.replace("GUT_GENOME", ""))
                        if gid_int % n == i:
                            selected_genomes.add(gid)
    except:
        tsprint(f"ERROR:  Genomes argument is not a list of genome ids or slices: {g}")
        raise
    return sorted(selected_genomes)


def import_uhgg_master(args):

    # Fetch table of contents from s3.
    # This will be read separately by each species build subcommand, so we make a local copy.
    local_toc = os.path.basename(outputs.genomes)
    command(f"rm -f {local_toc}")
    command(f"aws s3 cp --only-show-errors {outputs.genomes} {local_toc}")

    db = UHGG(local_toc)
    species_for_genome = db.genomes

    def genome_work(genome_id):
        assert genome_id in species_for_genome, f"Genome {genome_id} is not in the database."
        species_id = species_for_genome[genome_id]

        dest_file = imported_genome_file(genome_id, species_id, f"{genome_id}.fna.lz4")
        msg = f"Importing genome {genome_id} from species {species_id}."
        if find_files_with_retry(dest_file):
            if not args.force:
                tsprint(f"Destination {dest_file} for genome {genome_id} already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Importing", "Reimporting")

        tsprint(msg)
        slave_log = "import_uhgg.log"
        slave_subdir = f"{species_id}__{genome_id}"
        if not args.debug:
            command(f"rm -rf {slave_subdir}")
        if not os.path.isdir(slave_subdir):
            command(f"mkdir {slave_subdir}")
        # Recurisve call via subcommand.  Use subdir, redirect logs.
        slave_cmd = f"cd {slave_subdir}; PYTHONPATH={pythonpath()} {sys.executable} -m iggtools import_uhgg --genome {genome_id} --zzz_slave_mode --zzz_slave_toc {os.path.abspath(local_toc)} {'--debug' if args.debug else ''} &>> {slave_log}"
        with open(f"{slave_subdir}/{slave_log}", "w") as slog:
            slog.write(msg + "\n")
            slog.write(slave_cmd + "\n")
        try:
            command(slave_cmd)
        finally:
            # Cleanup should not raise exceptions of its own, so as not to interfere with any
            # prior exceptions that may be more informative.  Hence check=False.
            upload(f"{slave_subdir}/{slave_log}", imported_genome_file(genome_id, species_id, slave_log + ".lz4"), check=False)
            if not args.debug:
                command(f"rm -rf {slave_subdir}", check=False)

    genome_id_list = decode_genomes_arg(args, species_for_genome)
    multithreading_map(genome_work, genome_id_list, num_threads=CONCURRENT_GENOME_IMPORTS)


def import_uhgg_slave(args):
    """
    https://github.com/czbiohub/iggtools/wiki
    """

    violation = "Please do not call build_pangenome_slave directly.  Violation"
    assert args.zzz_slave_mode, f"{violation}:  Missing --zzz_slave_mode arg."
    assert os.path.isfile(args.zzz_slave_toc), f"{violation}: File does not exist: {args.zzz_slave_toc}"

    db = UHGG(args.zzz_slave_toc)
    representatives = db.representatives
    species_for_genome = db.genomes

    genome_id = args.genomes
    species_id = species_for_genome[genome_id]
    representative_id = representatives[species_id]

    dest = imported_genome_file(genome_id, species_id, f"{genome_id}.fna.lz4")
    command(f"aws s3 rm --recursive {imported_genome_file(genome_id, species_id, '')}")
    cleaned = clean_genome(genome_id, representative_id)
    upload(cleaned, dest)


def register_args(main_func):
    subparser = add_subcommand('import_uhgg', main_func, help='import selected genomes from UHGG')
    subparser.add_argument('--genomes',
                           dest='genomes',
                           required=False,
                           help="genome[,genome...] to import;  alternatively, slice in format idx:modulus, e.g. 1:30, meaning import genomes whose ids are 1 mod 30; or, the special keyword 'all' meaning all genomes")
    subparser.add_argument('--zzz_slave_toc',
                           dest='zzz_slave_toc',
                           required=False,
                           help=SUPPRESS) # "reserved to pass table of contents from master to slave"
    return main_func


@register_args
def main(args):
    tsprint(f"Executing iggtools subcommand {args.subcommand} with args {vars(args)}.")
    import_uhgg(args)
