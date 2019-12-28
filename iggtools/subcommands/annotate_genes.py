import os
import sys
from multiprocessing import Semaphore
from iggtools.common.argparser import add_subcommand, SUPPRESS
from iggtools.common.utils import tsprint, retry, command, multithreading_map, find_files, upload, pythonpath, upload_star, download_reference
from iggtools.models.uhgg import UHGG
from iggtools.params import outputs


CONCURRENT_PROKKA_RUNS = Semaphore(6)


def imported_genome_file(genome_id, species_id, component):
    return f"{outputs.cleaned_imports}/{species_id}/{genome_id}/{component}"


def annotations_file(genome_id, species_id, filename):
    # s3://microbiome-igg/2.0/prodigal/GUT_GENOMEDDDDDD.{fna, faa, gff, log}
    return f"{outputs.annotations}/{species_id}/{genome_id}/{filename}"


def unified_genome_id(genome_id):
    return "UHGG" + genome_id.replace("GUT_GENOME", "")


# 1. Occasional failures in aws s3 cp require a retry.
@retry
def download_genome(genome_id, cleaned_genome):
    command(f"rm -f {genome_id}.fasta")
    command(f"aws s3 cp --only-show-errors {cleaned_genome} - | lz4 -dc > {genome_id}.fasta")


def annotate_genome(genome_id, species_id):
    # Prokka will crash if installed <6 months ago.  It's a feature.  See tbl2asn.
    cleaned_genome = imported_genome_file(genome_id, species_id, f"{genome_id}.fna.lz4")
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


def annotate_genes(args):
    if args.zzz_slave_toc:
        annotate_genes_slave(args)
    else:
        annotate_genes_master(args)


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


def annotate_genes_master(args):

    # Fetch table of contents from s3.
    # This will be read separately by each species build subcommand, so we make a local copy.
    local_toc = download_reference(outputs.genomes)

    db = UHGG(local_toc)
    species_for_genome = db.genomes

    def genome_work(genome_id):
        assert genome_id in species_for_genome, f"Genome {genome_id} is not in the database."
        species_id = species_for_genome[genome_id]

        dest_file = annotations_file(genome_id, species_id, f"{genome_id}.fna.lz4")
        msg = f"Annotating genome {genome_id} from species {species_id}."
        if find_files_with_retry(dest_file):
            if not args.force:
                tsprint(f"Destination {dest_file} for genome {genome_id} annotations already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Importing", "Reimporting")


        with CONCURRENT_PROKKA_RUNS:

            tsprint(msg)
            slave_log = "annotate_genes.log"
            slave_subdir = f"{species_id}__{genome_id}"
            if not args.debug:
                command(f"rm -rf {slave_subdir}")
            if not os.path.isdir(slave_subdir):
                command(f"mkdir {slave_subdir}")
            # Recurisve call via subcommand.  Use subdir, redirect logs.
            slave_cmd = f"cd {slave_subdir}; PYTHONPATH={pythonpath()} {sys.executable} -m iggtools annotate_genes --genome {genome_id} --zzz_slave_mode --zzz_slave_toc {os.path.abspath(local_toc)} {'--debug' if args.debug else ''} &>> {slave_log}"
            with open(f"{slave_subdir}/{slave_log}", "w") as slog:
                slog.write(msg + "\n")
                slog.write(slave_cmd + "\n")
            try:
                command(slave_cmd)
            finally:
                # Cleanup should not raise exceptions of its own, so as not to interfere with any
                # prior exceptions that may be more informative.  Hence check=False.
                upload(f"{slave_subdir}/{slave_log}", annotations_file(genome_id, species_id, slave_log + ".lz4"), check=False)
                if not args.debug:
                    command(f"rm -rf {slave_subdir}", check=False)

    genome_id_list = decode_genomes_arg(args, species_for_genome)
    multithreading_map(genome_work, genome_id_list, num_threads=20)


def drop_lz4(filename):
    assert filename.endswith(".lz4")
    return filename[:-4]


def annotate_genes_slave(args):
    """
    https://github.com/czbiohub/iggtools/wiki
    """

    violation = "Please do not call build_pangenome_slave directly.  Violation"
    assert args.zzz_slave_mode, f"{violation}:  Missing --zzz_slave_mode arg."
    assert os.path.isfile(args.zzz_slave_toc), f"{violation}: File does not exist: {args.zzz_slave_toc}"

    db = UHGG(args.zzz_slave_toc)
    species_for_genome = db.genomes

    genome_id = args.genomes
    species_id = species_for_genome[genome_id]

    last_output = f"{genome_id}.fna.lz4"
    dest_file = annotations_file(genome_id, species_id, last_output)
    command(f"aws s3 rm --recursive {dest_file.rsplit('/', 1)[0]}")
    output_files = annotate_genome(genome_id, species_id)
    upload_tasks = []
    for o in output_files:
        olz = o + ".lz4"
        if olz != last_output:
            upload_tasks.append((o, annotations_file(genome_id, species_id, olz)))

    multithreading_map(upload_star, upload_tasks)

    # Upload this last because it indicates all other work has succeeded.
    upload(drop_lz4(last_output), annotations_file(genome_id, species_id, last_output))


def register_args(main_func):
    subparser = add_subcommand('annotate_genes', main_func, help='annotate selected genomes with prokka')
    subparser.add_argument('--genomes',
                           dest='genomes',
                           required=False,
                           help="genome[,genome...] to import;  alternatively, slice in format idx:modulus, e.g. 1:30, meaning annotate genomes whose ids are 1 mod 30; or, the special keyword 'all' meaning all genomes")
    subparser.add_argument('--zzz_slave_toc',
                           dest='zzz_slave_toc',
                           required=False,
                           help=SUPPRESS) # "reserved to pass table of contents from master to slave"
    return main_func


@register_args
def main(args):
    tsprint(f"Executing iggtools subcommand {args.subcommand} with args {vars(args)}.")
    annotate_genes(args)
