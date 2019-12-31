import json
from collections import defaultdict
import numpy as np
import Bio.SeqIO
import pysam

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, num_physical_cores, command, InputStream, OutputStream, select_from_tsv, multithreading_hashmap, download_reference, split
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
                           help='Alignment speed/sensitivity (very-sensitive).  If aln_mode is local (default) this automatically issues the corresponding very-sensitive-local, etc flag to bowtie2.')
    subparser.add_argument('--aln_mode',
                           type=str,
                           dest='aln_mode',
                           default='local',
                           choices=['local', 'global'],
                           help='Global/local read alignment (local, corresponds to the bowtie2 --local).  Global corresponds to the bowtie2 default --end-to-end.')
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
    with InputStream(f"{outdir}/species/species_profile.txt") as stream:
        return dict(select_from_tsv(stream, {"species_id": str, "coverage": float}))


def select_species(species_profile, coverage_threshold):
    return {species_id: species_coverage for species_id, species_coverage in species_profile.items() if species_coverage >= coverage_threshold}


def pangenome_file(species_id, component):
    # s3://microbiome-igg/2.0/pangenomes/GUT_GENOMEDDDDDD/{genes.ffn, centroids.ffn, gene_info.txt}
    return f"{outputs.pangenomes}/{species_id}/{component}"


def pangenome_align(args, tempdir):
    """ Use Bowtie2 to map reads to all specified genome species """

    max_reads = f"-u {args.max_reads}" if args.max_reads else ""
    aln_mode = "local" if args.aln_mode == "local" else "end-to-end"
    aln_speed = args.aln_speed if aln_mode == "end_to_end" else args.aln_speed + "-local"
    r2 = ""
    if args.r2:
        r1 = f"-1 {args.r1}"
        r2 = f"-2 {args.r2}"
    elif args.aln_interleaved:
        r1 = f"--interleaved {args.r1}"
    else:
        r1 = f"-U {args.r1}"

    try:
        command(f"set -o pipefail;  bowtie2 --no-unal -x {tempdir}/pangenomes {max_reads} --{aln_mode} --{aln_speed} --threads {num_physical_cores} -q {r1} {r2} | samtools view --threads {num_physical_cores} -b - > {tempdir}/pangenomes.bam")
    except:
        command(f"rm -f {tempdir}/pangenomes.bam")
        raise


def keep_read(aln, min_pid, min_readq, min_mapq, min_aln_cov):
    align_len = len(aln.query_alignment_sequence)
    query_len = aln.query_length
    # min pid
    if 100 * (align_len - dict(aln.tags)['NM']) / float(align_len) < min_pid:
        return False
    # min read quality
    if np.mean(aln.query_qualities) < min_readq:
        return False
    # min map quality
    if aln.mapping_quality < min_mapq:
        return False
    # min aln cov
    if align_len / float(query_len) < min_aln_cov:
        return False
    return True


def count_mapped_bp(args, tempdir, genes):
    """ Count number of bp mapped to each gene across pangenomes """
    bam_path = f"{tempdir}/genes/temp/pangenomes.bam'"
    bamfile = pysam.AlignmentFile(bam_path, "rb")

    # loop over alignments, sum values per gene
    for aln in bamfile.fetch(until_eof=True):
        gene_id = bamfile.getrname(aln.reference_id)
        gene = genes[gene_id]
        gene["aligned_reads"] += 1
        if keep_read(aln, args.mapid, args.readq, args.mapq, args.aln_cov):
            gene["mapped_reads"] += 1
            gene["depth"] += len(aln.query_alignment_sequence) / float(gene["length"])

    tsprint("Pangenome count_mapped_bp:  total aligned reads: %s" % sum(g["aligned_reads"] for g in genes.values()))
    tsprint("Pangenome count_mapped_bp:  total mapped reads: %s" % sum(g["mapped_reads"] for g in genes.values()))

    # Group gene depths by species
    gene_depths = defaultdict(list)
    for g in genes.values():
        gene_depths[g["species_id"]].append(g["depth"])

    # loop over species, compute summaries
    num_covered_genes = {}
    mean_coverage = {}
    for species_id, depths in gene_depths.items():
        non_zero_depths = [d for d in depths if d > 0]
        num_covered_genes = len(non_zero_depths)
        mean_coverage[species_id] = np.mean(non_zero_depths) if non_zero_depths else 0
        num_covered_genes[species_id] = num_covered_genes

    return num_covered_genes, mean_coverage


def normalize(genes):
    """ Count number of bp mapped to each marker gene """
    # compute marker depth
    species_markers = defaultdict(lambda: defaultdict(int))
    for gene in genes.values():
        gene_marker_id = gene["marker_id"]
        if gene_marker_id:
            species_markers[gene["species_id"]][gene_marker_id] += gene["depth"]
    # compute median marker depth
    species_markers_coverage = {species_id: np.median(marker_depths) for species_id, marker_depths in species_markers.items()}
    # normalize genes by median marker depth
    for gene in genes.values():
        species_id = gene["species_id"]
        marker_coverage = species_markers_coverage[species_id]
        if marker_coverage > 0:
            gene["copies"] = gene["depth"] / marker_coverage
    return species_markers_coverage


def write_results(tempdir, species, num_covered_genes, markers_coverage, mean_coverage):
    # open outfiles for each species_id
    header = ['gene_id', 'count_reads', 'coverage', 'copy_number']
    for species_id, species_genes in species:
        path = f"{tempdir}/genes/output/{species_id}.genes.gz"
        with OutputStream(path) as sp_out:
            sp_out.write('\t'.join(header) + '\n')
            for gene_id, gene in species_genes.items():
                values = [gene_id, gene["mapped_reads"], gene["depth"], gene["copies"]]
                sp_out.write('\t'.join([str(v) for v in values]) + '\n')
    # summary stats
    header = ['species_id', 'pangenome_size', 'covered_genes', 'fraction_covered', 'mean_coverage', 'marker_coverage', 'aligned_reads', 'mapped_reads']
    path = f"{tempdir}/genes/summary.txt"
    with OutputStream(path) as file:
        file.write('\t'.join(header) + '\n')
        for species_id, species_genes in species.items():
            aligned_reads = sum(g["algined_reads"] for g in species_genes)
            mapped_reads = sum(g["mapped_reads"] for g in species_genes)
            pangenome_size = len(species_genes)
            values = [species_id, pangenome_size, num_covered_genes[species_id], num_covered_genes[species_id] / pangenome_size, mean_coverage[species_id], markers_coverage[species_id], aligned_reads, mapped_reads]
            file.write('\t'.join(str(v) for v in values) + '\n')


def pangenome_coverage(args, tempdir, species, genes):
    """ Compute coverage of pangenome for species_id and write results to disk """
    num_covered_genes, species_mean_coverage = count_mapped_bp(args, tempdir, genes)
    species_markers_coverage = normalize(genes)
    write_results(tempdir, species, num_covered_genes, species_markers_coverage, species_mean_coverage)


def scan_centroids(centroids_files):
    species = defaultdict(dict)
    genes = {}
    for species_id, centroid_filename in centroids_files:
        with InputStream(centroid_filename) as file:
            for centroid in Bio.SeqIO.parse(file, 'fasta'):
                centroid_gene_id = centroid.id
                centroid_gene = {
                    "centroid_gene_id": centroid_gene_id,
                    "species_id": species_id,
                    "length": len(centroid.seq),
                    "depth": 0.0,
                    "aligned_reads": 0,
                    "mapped_reads": 0,
                    "copies": 0.0,
                }
                species[species_id][centroid_gene_id] = centroid_gene
                genes[centroid_gene_id] = centroid_gene
    return species, genes


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
    centroids_files = multithreading_hashmap(download_centroid, species_profile.keys(), num_threads=20)

    build_pangenome_db(tempdir, centroids_files)

    pangenome_align(args, tempdir)

    species, genes = scan_centroids(centroids_files)

    pangenome_coverage(args, tempdir, species, genes)


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_run_genes(args)
