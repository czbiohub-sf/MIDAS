import json
from collections import defaultdict

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, num_physical_cores, command, InputStream, OutputStream, parse_table, multithreading_hashmap, download_reference, split
from iggtools.params import inputs, outputs


DEFAULT_ALN_COV = 0.75
DEFAULT_SPECIES_COVERAGE = 3.0


def register_args(main_func):
    subparser = add_subcommand('midas_run_genes', main_func, help='metagenomic pan-genome profiling')
    subparser.add_argument('outdir',
                           type=str,
                           help="""Path to directory to store results.  Name should correspond to unique sample identifier.""")
    subparser.add_argument('-1',
                           dest='r1',
                           required=True,
                           help="FASTA/FASTQ file containing 1st mate if using paired-end reads.  Otherwise FASTA/FASTQ containing unpaired reads.")
    subparser.add_argument('-2',
                           dest='r2',
                           help="FASTA/FASTQ file containing 2nd mate if using paired-end reads.")
    subparser.add_argument('--mapid',
                           dest='mapid',
                           type=float,
                           metavar="FLOAT",
                           help=f"Discard reads with alignment identity < MAPID.  Values between 0-100 accepted.  By default gene-specific species-level cutoffs are used, as specifeid in {inputs.marker_genes_hmm_cutoffs}")
    subparser.add_argument('--max_reads',
                           dest='max_reads',
                           type=int,
                           metavar="INT",
                           help=f"Number of reads to use from input file(s).  (All)")
    subparser.add_argument('--species_cov',
                           type=float,
                           dest='species_cov',
                           metavar='FLOAT',
                           default=DEFAULT_SPECIES_COVERAGE,
                           help=f"Include species with >X coverage ({DEFAULT_SPECIES_COVERAGE})")
    if False:
        # This is unused.
        subparser.add_argument('--species_topn',
                               type=int,
                               dest='species_topn',
                               metavar='INT',
                               help='Include top N most abundant species')

    #  Alignment flags (bowtie, or postprocessing)
    subparser.add_argument('--aln_cov',
                           dest='aln_cov',
                           default=DEFAULT_ALN_COV,
                           type=float,
                           metavar="FLOAT",
                           help=f"Discard reads with alignment coverage < ALN_COV ({DEFAULT_ALN_COV}).  Values between 0-1 accepted.")
    subparser.add_argument('--readq',
                           dest='readq',
                           type=int,
                           metavar="INT",
                           help=f"Discard reads with mean quality < READQ (20)")
    subparser.add_argument('--aln_speed',
                           type=str,
                           dest='aln_speed',
                           default='very-sensitive',
                           choices=['very-fast', 'fast', 'sensitive', 'very-sensitive'],
                           help='Alignment speed/sensitivity (very-sensitive)')
    subparser.add_argument('--aln_mode',
                           type=str,
                           dest='aln_mode',
                           default='local',
                           choices=['local', 'global'],
                           help='Global/local read alignment (local)')
    subparser.add_argument('--aln_interleaved',
                           action='store_true',
                           default=False,
                           help='FASTA/FASTQ file in -1 are paired and contain forward AND reverse reads')

    return main_func


def build_pangenome_db(tempdir, centroids):
    command(f"touch {tempdir}/pangenomes.fa")
    for files in split(centroids.values(), 20):  # keep "cat" commands short
        command("cat " + " ".join(files) + f" >> {tempdir}/pangenomes.fa")
    command(f"bowtie2-build --threads {num_physical_cores} {tempdir}/pangenomes.fa {tempdir}/pangenomes")


def parse_species_profile(outdir):
    "Return map of species_id to coverage for the species present in the sample."
    spfilename = f"{outdir}/species/species_profile.txt"
    with InputStream(spfilename) as sppf:
        return {species_id: float(species_coverage_str) for species_id, species_coverage_str in parse_table(sppf, ["species_id", "coverage"])}


def select_species(species_profile, coverage_threshold):
    return {species_id: species_coverage for species_id, species_coverage in species_profile.items() if species_coverage > coverage_threshold}


def pangenome_file(species_id, component):
    # s3://microbiome-igg/2.0/pangenomes/GUT_GENOMEDDDDDD/{genes.ffn, centroids.ffn, gene_info.txt}
    return f"{outputs.pangenomes}/{species_id}/{component}"


def pangenome_align(args, tempdir):
    """ Use Bowtie2 to map reads to all specified genome species """

    max_reads = f"-u {args.max_reads}" if args.max_reads else ""
    local_vs_global = "--local" if args.aln_mode == "local" else ""
    r2 = ""
    if args.r2:
        r1 = f"-1 {args.r1}"
        r2 = f"-2 {args.r2}"
    elif args.aln_interleaved:
        r1 = f"--interleaved {args.r1}"
    else:
        r1 = f"-U {args.r1}"

    try:
        command(f"set -o pipefail;  bowtie2 --no-unal -x {tempdir}/pangenomes {max_reads} --{args.aln_speed} {local_vs_global} --threads {num_physical_cores} -q {r1} {r2} | samtools view --threads {num_physical_cores} -b - > {tempdir}/pangenomes.bam")
    except:
        command(f"rm -f {tempdir}/pangenomes.bam")
        raise


def midas_run_genes(args):

    tempdir = f"{args.outdir}/genes/temp_sc{args.species_cov}"

    command(f"rm -rf {tempdir}")
    command(f"mkdir -p {tempdir}")

    # Species profile must exist, it's output by run_midas_species
    full_species_profile = parse_species_profile(args.outdir)

    species_profile = select_species(full_species_profile, args.species_cov)

    def download_centroid(species_id):
        return download_reference(pangenome_file(species_id, "centroids.ffn.lz4"), f"{tempdir}/{species_id}")  # TODO colocate samples to overlap reference downloads

    # Download centroids.ffn for every species
    centroids = multithreading_hashmap(download_centroid, species_profile.keys(), num_threads=20)

    build_pangenome_db(tempdir, centroids)

    pangenome_align(args, tempdir)


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_run_genes(args)
