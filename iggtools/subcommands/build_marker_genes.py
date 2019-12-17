import os
import sys
from collections import defaultdict
import Bio.SeqIO
from iggtools.common.argparser import add_subcommand, SUPPRESS
from iggtools.common.utils import tsprint, InputStream, parse_table, retry, command, multithreading_map, find_files, sorted_dict, upload, pythonpath
from iggtools.params import inputs, outputs


CONCURRENT_MARKER_GENES_IDENTIFY = 20


def input_annotations_file(genome_id, species_id, filename):
    # s3://microbiome-igg/2.0/gene_annotations/{SPECIES_ID}/{GENOME_ID}/{GENOME_ID}.{fna, faa, gff, log}
    return f"{outputs.annotations}/{species_id}/{genome_id}/{filename}"


def output_marker_genes_file(genome_id, species_id, filename):
    # s3://{igg}/marker_genes/phyeco/temp/GUT_GENOME138501.{hmmsearch, markers.fa, markers.map}
    return f"{outputs.marker_genes}/temp/{species_id}/{genome_id}/{filename}"


def read_toc(genomes_tsv, deep_sort=False):
    # Read in the table of contents.
    # We will centralize and improve this soon.
    species = defaultdict(dict)
    representatives = {}
    genomes = {}
    with InputStream(genomes_tsv) as table_of_contents:
        for row in parse_table(table_of_contents, ["genome", "species", "representative", "genome_is_representative"]):
            genome_id, species_id, representative_id, _ = row
            species[species_id][genome_id] = row
            representatives[species_id] = representative_id
            genomes[genome_id] = species_id
    if deep_sort:
        for sid in species.keys():
            species[sid] = sorted_dict(species[sid])
        species = sorted_dict(species)
    return species, representatives, genomes


# 1. Occasional failures in aws s3 cp require a retry.
@retry
def download_genome(genome_id, annotated_genes):
    command(f"rm -f {genome_id}.faa")
    command(f"aws s3 cp --only-show-errors {annotated_genes} - | lz4 -dc > {genome_id}.faa")

@retry
def download_marker_genes_hmm():
    command(f"rm -f marker_genes.hmm")
    command(f"aws s3 cp --only-show-errors {inputs.marker_genes_hmm} marker_genes.hmm")


def hmmsearch(genome_id, species_id, num_threads=1):
    # Input
    annotated_genes = input_annotations_file(genome_id, species_id, f"{genome_id}.faa.lz4")
    download_genome(genome_id, annotated_genes)
    download_marker_genes_hmm()

    # Output
    marker_hmmsearch = f"{genome_id}.hmmsearch"

    # Command
    if find_files(marker_hmmsearch):
        tsprint(f"Found hmmsearch results for genome {genome_id} from prior run.")
    else:
        try:
            command(f"hmmsearch --noali --cpu {num_threads} --domtblout {marker_hmmsearch} marker_genes.hmm {genome_id}.faa")
        except:
            # Do not keep bogus zero-length files;  those are harmful if we rerun in place.
            command(f"mv {marker_hmmsearch} {marker_hmmsearch}.bogus", check=False)
            raise

    return marker_hmmsearch


def build_marker_genes(args):
    if args.zzz_slave_toc:
        build_marker_genes_slave(args)
    else:
        build_marker_genes_master(args)


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


def build_marker_genes_master(args):

    # Fetch table of contents from s3.
    # This will be read separately by each species build subcommand, so we make a local copy.
    local_toc = os.path.basename(outputs.genomes)
    command(f"rm -f {local_toc}")
    command(f"aws s3 cp --only-show-errors {outputs.genomes} {local_toc}")

    _, _, species_for_genome = read_toc(local_toc)

    def genome_work(genome_id):
        assert genome_id in species_for_genome, f"Genome {genome_id} is not in the database."
        species_id = species_for_genome[genome_id]

        dest_file = output_marker_genes_file(genome_id, species_id, f"{genome_id}.hmmsearch.lz4")
        msg = f"Running HMMsearch for genome {genome_id} from species {species_id}."
        if find_files_with_retry(dest_file):
            if not args.force:
                tsprint(f"Destination {dest_file} for genome {genome_id} already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Running", "Rerunning")

        tsprint(msg)
        slave_log = "build_marker_genes.log"
        slave_subdir = f"{species_id}__{genome_id}"
        if not args.debug:
            command(f"rm -rf {slave_subdir}")
        if not os.path.isdir(slave_subdir):
            command(f"mkdir {slave_subdir}")

        # Recurisve call via subcommand.  Use subdir, redirect logs.
        slave_cmd = f"cd {slave_subdir}; PYTHONPATH={pythonpath()} {sys.executable} -m iggtools build_marker_genes --genome {genome_id} --zzz_slave_mode --zzz_slave_toc {os.path.abspath(local_toc)} {'--debug' if args.debug else ''} &>> {slave_log}"
        with open(f"{slave_subdir}/{slave_log}", "w") as slog:
            slog.write(msg + "\n")
            slog.write(slave_cmd + "\n")
        try:
            command(slave_cmd)
        finally:
            # Cleanup should not raise exceptions of its own, so as not to interfere with any
            # prior exceptions that may be more informative.  Hence check=False.
            upload(f"{slave_subdir}/{slave_log}", output_marker_genes_file(genome_id, species_id, slave_log + ".lz4"), check=False)
            if not args.debug:
                command(f"rm -rf {slave_subdir}", check=False)

    genome_id_list = decode_genomes_arg(args, species_for_genome)
    multithreading_map(genome_work, genome_id_list, num_threads=CONCURRENT_MARKER_GENES_IDENTIFY)


def build_marker_genes_slave(args):
    """
    https://github.com/czbiohub/iggtools/wiki
    """

    violation = "Please do not call build_merker_genes_slave directly.  Violation"
    assert args.zzz_slave_mode, f"{violation}:  Missing --zzz_slave_mode arg."
    assert os.path.isfile(args.zzz_slave_toc), f"{violation}: File does not exist: {args.zzz_slave_toc}"

    _, _, species_for_genome = read_toc(args.zzz_slave_toc)

    genome_id = args.genomes
    species_id = species_for_genome[genome_id]

    dest_file = output_marker_genes_file(genome_id, species_id, f"{genome_id}.hmmsearch.lz4")
    command(f"aws s3 rm --recursive {output_marker_genes_file(genome_id, species_id, '')}")
    marker_hmmsearch = hmmsearch(genome_id, species_id, num_threads=1)
    upload(marker_hmmsearch, dest_file)


def register_args(main_func):
    subparser = add_subcommand('build_marker_genes', main_func, help='identify marker genes for  given genomes')
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
    build_marker_genes(args)
