import json
import os
from collections import defaultdict
import numpy as np
import Bio.SeqIO
from pysam import AlignmentFile  # pylint: disable=no-name-in-module

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, num_physical_cores, command, InputStream, OutputStream, select_from_tsv, multithreading_hashmap, download_reference
from iggtools.params import outputs
from iggtools.common.samples import parse_species_profile, select_species
from iggtools.common.bowtie2 import build_bowtie2_db


DEFAULT_ALN_COV = 0.75
DEFAULT_SPECIES_COVERAGE = 3.0
DEFAULT_ALN_MAPID = 94.0
DEFAULT_ALN_READQ = 20
DEFAULT_ALN_MAPQ = 0

DECIMALS = ".6f"

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
    subparser.add_argument('--aln_readq',
                           dest='aln_readq',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_ALN_READQ,
                           help=f"Discard reads with mean quality < READQ ({DEFAULT_ALN_READQ})")
    subparser.add_argument('--aln_mapid',
                           dest='aln_mapid',
                           type=float,
                           metavar="FLOAT",
                           default=DEFAULT_ALN_MAPID,
                           help=f"Discard reads with alignment identity < MAPID.  Values between 0-100 accepted.  ({DEFAULT_ALN_MAPID})")
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


def pangenome_file(species_id, component):
    # s3://microbiome-igg/2.0/pangenomes/GUT_GENOMEDDDDDD/{genes.ffn, centroids.ffn, gene_info.txt}
    return f"{outputs.pangenomes}/{species_id}/{component}"


def pangenome_align(args, tempdir):
    """ Use Bowtie2 to map reads to all specified genome species """

    if args.debug and os.path.exists(f"{tempdir}/pangenomes.bam"):
        tsprint(f"Skipping alignment in debug mode as temporary data exists: {tempdir}/pangenomes.bam")
        return

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
        command(f"set -o pipefail; bowtie2 --no-unal -x {tempdir}/pangenomes {max_reads} --{aln_mode} --{aln_speed} --threads {num_physical_cores} -q {r1} {r2} | samtools view --threads {num_physical_cores} -b - > {tempdir}/pangenomes.bam")
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
    """ Count number of bp mapped to each gene across pangenomes.  Return number covered genes and average gene depth per species.  Result contains only covered species, but being a defaultdict, would yield 0 for any uncovered species, which is appropriate. """
    bam_path = f"{tempdir}/pangenomes.bam"
    bamfile = AlignmentFile(bam_path, "rb")
    covered_genes = {}

    # loop over alignments, sum values per gene
    for aln in bamfile.fetch(until_eof=True):
        gene_id = bamfile.getrname(aln.reference_id)
        gene = genes[gene_id]
        gene["aligned_reads"] += 1
        if keep_read(aln, args.aln_mapid, args.aln_readq, DEFAULT_ALN_MAPQ, args.aln_cov):
            gene["mapped_reads"] += 1
            gene["depth"] += len(aln.query_alignment_sequence) / float(gene["length"])
            covered_genes[gene_id] = gene

    tsprint("Pangenome count_mapped_bp:  total aligned reads: %s" % sum(g["aligned_reads"] for g in genes.values()))
    tsprint("Pangenome count_mapped_bp:  total mapped reads: %s" % sum(g["mapped_reads"] for g in genes.values()))

    # Filter to genes with non-zero depth, then group by species
    nonzero_gene_depths = defaultdict(list)
    for g in covered_genes.values():
        gene_depth = g["depth"]
        if gene_depth > 0:  # This should always pass, because ags.aln_cov is always >0.
            species_id = g["species_id"]
            nonzero_gene_depths[species_id].append(gene_depth)

    # Compute number of covered genes per species, and average gene depth.
    num_covered_genes = defaultdict(int)
    mean_coverage = defaultdict(float)
    for species_id, non_zero_depths in nonzero_gene_depths.items():
        num_covered_genes[species_id] = len(non_zero_depths)
        mean_coverage[species_id] = np.mean(non_zero_depths)

    return num_covered_genes, mean_coverage, covered_genes


def normalize(genes, covered_genes, markers):
    """ Normalize gene depth by median marker depth, to infer gene copy count.  Return median marker coverage for each covered species in a defaultdict, so that accessing an uncovered species would appropriately yield 0. """
    # compute marker depth
    species_markers = defaultdict(lambda: defaultdict(int))
    for gene_id, marker_id in markers:
        g = genes[gene_id]
        species_markers[g["species_id"]][marker_id] += g["depth"]
    # compute median marker depth for each species
    species_markers_coverage = defaultdict(float)
    for species_id, marker_depths in species_markers.items():
        species_markers_coverage[species_id] = np.median(list(marker_depths.values()))   # np.median doesn't take iterators
    # infer copy count for each covered gene
    for gene in covered_genes.values():
        species_id = gene["species_id"]
        marker_coverage = species_markers_coverage[species_id]
        if marker_coverage > 0:
            gene["copies"] = gene["depth"] / marker_coverage
    return species_markers_coverage


def write_results(outdir, species, num_covered_genes, species_markers_coverage, species_mean_coverage):
    if not os.path.exists(f"{outdir}/genes/output"):
        command(f"mkdir -p {outdir}/genes/output")
    # open outfiles for each species_id
    header = ['gene_id', 'count_reads', 'coverage', 'copy_number']
    for species_id, species_genes in species.items():
        path = f"{outdir}/genes/output/{species_id}.genes.lz4"
        with OutputStream(path) as sp_out:
            sp_out.write('\t'.join(header) + '\n')
            for gene_id, gene in species_genes.items():
                if gene["depth"] == 0:
                    # Sparse by default here.  You can get the pangenome_size from the summary file, emitted below.
                    continue
                values = [gene_id, str(gene["mapped_reads"]), format(gene["depth"], DECIMALS), format(gene["copies"], DECIMALS)]
                sp_out.write('\t'.join(values) + '\n')
    # summary stats
    header = ['species_id', 'pangenome_size', 'covered_genes', 'fraction_covered', 'mean_coverage', 'marker_coverage', 'aligned_reads', 'mapped_reads']
    path = f"{outdir}/genes/summary.txt"
    with OutputStream(path) as file:
        file.write('\t'.join(header) + '\n')
        for species_id, species_genes in species.items():
            # No sparsity here -- should be extremely rare for a species row to be all 0.
            aligned_reads = sum(g["aligned_reads"] for g in species_genes.values())
            mapped_reads = sum(g["mapped_reads"] for g in species_genes.values())
            pangenome_size = len(species_genes)
            values = [
                species_id,
                str(pangenome_size),
                str(num_covered_genes[species_id]),
                format(num_covered_genes[species_id] / pangenome_size, DECIMALS),
                format(species_mean_coverage[species_id], DECIMALS),
                format(species_markers_coverage[species_id], DECIMALS),
                str(aligned_reads),
                str(mapped_reads)
            ]
            file.write('\t'.join(values) + '\n')


def scan_centroids(centroids_files):
    species = defaultdict(dict)
    genes = {}
    for species_id, centroid_filename in centroids_files.items():
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


def scan_markers(genes, marker_genes_map_file):
    markers = []
    with InputStream(marker_genes_map_file) as mg_map:
        for gene_id, marker_id in select_from_tsv(mg_map, ["gene_id", "marker_id"], {"species_id": str, "genome_id": str, "gene_id": str, "gene_len": int, "marker_id": str}):
            if gene_id in genes:
                markers.append((gene_id, marker_id))
    return markers


def midas_run_genes(args):

    tempdir = f"{args.outdir}/genes/temp_sc{args.species_cov}"

    if args.debug and os.path.exists(tempdir):
        tsprint(f"INFO:  Reusing existing temp data in {tempdir} according to --debug flag.")
    else:
        command(f"rm -rf {tempdir}")
        command(f"mkdir -p {tempdir}")

    try:
        # The full species profile must exist -- it is output by run_midas_species.
        # Restrict to species above requested coverage.
        full_species_profile = parse_species_profile(args.outdir)
        species_profile = select_species(full_species_profile, args.species_cov)

        def download_centroid(species_id):
            return download_reference(pangenome_file(species_id, "centroids.ffn.lz4"), f"{tempdir}/{species_id}")  # TODO colocate samples to overlap reference downloads

        # Download centroids.ffn for every species in the restricted species profile.
        centroids_files = multithreading_hashmap(download_centroid, species_profile.keys(), num_threads=20)

        # Perhaps avoid this giant conglomerated file, fetching instead submaps for each species.
        # Also colocate/cache/download in master for multiple slave subcommand invocations.
        marker_genes_map = "s3://microbiome-igg/2.0/marker_genes/phyeco/phyeco.map.lz4"
        build_bowtie2_db(tempdir, "pangenomes", centroids_files)
        pangenome_align(args, tempdir)

        # Compute coverage of pangenome for each present species and write results to disk
        species, genes = scan_centroids(centroids_files)
        num_covered_genes, species_mean_coverage, covered_genes = count_mapped_bp(args, tempdir, genes)
        markers = scan_markers(genes, marker_genes_map)
        species_markers_coverage = normalize(genes, covered_genes, markers)
        write_results(args.outdir, species, num_covered_genes, species_markers_coverage, species_mean_coverage)
    except:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error.  Specify --debug flag to keep.")
            command(f"rm -rf {tempdir}", check=False)


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_run_genes(args)
