import json
import os
from collections import defaultdict
import numpy as np
import Bio.SeqIO
from pysam import AlignmentFile  # pylint: disable=no-name-in-module

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, command, InputStream, OutputStream, select_from_tsv, multithreading_hashmap, download_reference
from iggtools.params import outputs
from iggtools.common.bowtie2 import build_bowtie2_db, bowtie2_align
from iggtools.models.uhgg import UHGG, pangenome_file
from iggtools.params.outputs import marker_genes
from iggtools.params.schemas import genes_profile_schema, genes_info_schema, genes_schema, MARKER_INFO_SCHEMA


DEFAULT_ALN_COV = 0.75
DEFAULT_GENOME_COVERAGE = 3.0
DEFAULT_ALN_MAPID = 94.0
DEFAULT_ALN_READQ = 20
DEFAULT_ALN_MAPQ = 0

DECIMALS = ".6f"

def register_args(main_func):
    subparser = add_subcommand('midas_run_genes', main_func, help='metagenomic pan-genome profiling')
    subparser.add_argument('midas_outdir',
                           type=str,
                           help="""Path to directory to store results.  Name should correspond to unique sample identifier.""")
    subparser.add_argument('-1',
                           dest='r1',
                           required=True,
                           help="FASTA/FASTQ file containing 1st mate if using paired-end reads.  Otherwise FASTA/FASTQ containing unpaired reads.")
    subparser.add_argument('-2',
                           dest='r2',
                           help="FASTA/FASTQ file containing 2nd mate if using paired-end reads.")

    subparser.add_argument('--genome_coverage',
                           type=float,
                           dest='genome_coverage',
                           metavar='FLOAT',
                           default=DEFAULT_GENOME_COVERAGE,
                           help=f"Include species with >X coverage ({DEFAULT_GENOME_COVERAGE})")
    subparser.add_argument('--species_list',
                           dest='species_list',
                           type=str,
                           metavar="CHAR",
                           help=f"Comma separated list of species ids")
    subparser.add_argument('--bt2_db_indexes',
                           dest='bt2_db_indexes',
                           type=str,
                           metavar="CHAR",
                           help=f"Prebuilt bowtie2 database indexes")

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
    subparser.add_argument('--aln_mapq',
                           dest='aln_mapq',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_ALN_MAPQ,
                           help=f"Discard reads with DEFAULT_ALN_MAPQ < MAPQ. ({DEFAULT_ALN_MAPQ})")
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

    subparser.add_argument('--max_reads',
                           dest='max_reads',
                           type=int,
                           metavar="INT",
                           help=f"Number of reads to use from input file(s).  (All)")
    return main_func


def format_data(x):
    return format(x, DECIMALS) if isinstance(x, float) else str(x)


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


def count_mapped_bp(pangenome_bamfile, genes):
    """ Count number of bp mapped to each gene across pangenomes.
    Return number covered genes and average gene depth per species.
    Result contains only covered species, but being a defaultdict,
    would yield 0 for any uncovered species, which is appropriate.
    """

    global global_args
    args = global_args

    #pangenome_bamfile = f"{tempdir}/pangenomes.bam" #<- pass this as an input
    bamfile = AlignmentFile(pangenome_bamfile, "rb")
    covered_genes = {}

    # loop over alignments, sum values per gene
    for aln in bamfile.fetch(until_eof=True):
        gene_id = bamfile.getrname(aln.reference_id)
        gene = genes[gene_id]
        gene["aligned_reads"] += 1
        if keep_read(aln, args.aln_mapid, args.aln_readq, args.aln_mapq, args.aln_cov):
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
    """ Normalize gene depth by median marker depth, to infer gene copy count.
    Return median marker coverage for each covered species in a defaultdict,
    so that accessing an uncovered species would appropriately yield 0.
    """
    # compute marker depth
    species_markers = defaultdict(lambda: defaultdict(int))
    for gene_id, marker_id in markers:
        g = genes[gene_id]
        species_markers[g["species_id"]][marker_id] += g["depth"]
    # compute median marker depth for each species
    species_markers_coverage = defaultdict(float)
    for species_id, marker_depths in species_markers.items():
        species_markers_coverage[species_id] = np.median(list(marker_depths.values())) # np.median doesn't take iterators
    # infer copy count for each covered gene
    for gene in covered_genes.values():
        species_id = gene["species_id"]
        marker_coverage = species_markers_coverage[species_id]
        if marker_coverage > 0:
            gene["copies"] = gene["depth"] / marker_coverage
    return species_markers_coverage


def write_results(outdir, species, num_covered_genes, species_markers_coverage, species_mean_coverage):
    # We should not handle this here
    #if not os.path.exists(f"{outdir}/genes/output"):
    #    command(f"mkdir -p {outdir}/genes/output")

    # open outfiles for each species_id
    header = list(genes_schema.keys())
    for species_id, species_genes in species.items():
        path = sample.layout(species_id)["genes_coverage"]
        #path = f"{outdir}/genes/output/{species_id}.genes.lz4"
        with OutputStream(path) as sp_out:
            sp_out.write('\t'.join(header) + '\n')
            for gene_id, gene in species_genes.items():
                if gene["depth"] == 0:
                    # Sparse by default here.  You can get the pangenome_size from the summary file, emitted below.
                    continue
                values = [gene_id, str(gene["mapped_reads"]), format(gene["depth"], DECIMALS), format(gene["copies"], DECIMALS)]
                sp_out.write('\t'.join(values) + '\n')
    # summary stats
    header = list(genes_profile_schema.keys())
    #header = ['species_id', 'pangenome_size', 'covered_genes', 'fraction_covered', 'mean_coverage', 'marker_coverage', 'aligned_reads', 'mapped_reads']
    #path = f"{outdir}/genes/summary.txt"
    path = sample.layout()["genes_summary"]
    with OutputStream(path) as file:
        file.write('\t'.join(header) + '\n')
        for species_id, species_genes in species.items():
            # No sparsity here -- should be extremely rare for a species row to be all 0.
            aligned_reads = sum(g["aligned_reads"] for g in species_genes.values())
            mapped_reads = sum(g["mapped_reads"] for g in species_genes.values())
            pangenome_size = len(species_genes)

            values = [
                species_id,
                pangenome_size,
                num_covered_genes[species_id],
                num_covered_genes[species_id] / pangenome_size,
                species_mean_coverage[species_id],
                species_markers_coverage[species_id],
                aligned_reads,
                mapped_reads
            ]

            file.write("\t".join(map(format_data, values)) + "\n")



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
        for gene_id, marker_id in select_from_tsv(mg_map, ["gene_id", "marker_id"], schema = MARKER_INFO_SCHEMA):
            if gene_id in genes:
                markers.append((gene_id, marker_id))
    return markers


def midas_run_genes(args):

    sample = Sample(args.sample_name, args.midas_outdir, "snps")
    sample.create_output_dir(args.debug)

    global global_sample
    global_sample = sample
    global global_args
    global_args = args

    try:
        if args.bt2_db_indexes:
            if bowtie2_index_exists(os.path.dirname(args.bt2_db_indexes), os.path.basename(bt2_db_indexes)):
                print("good the pre built bt2 index exist")
                bt2_db_dir = os.path.dirname(args.bt2_db_indexes)
                bt2_db_name = os.path.basename(bt2_db_indexes)
            else:
                print("Error: good the pangenome pre built bt2 index exist")
                exit(0)
        else:
            # add this to the mapping layout
            sample.bt2_db_dir = f"{sample.tempdir}/dbs"
            command(f"mkdir -p {sample.bt2_db_dir}")

            bt2_db_dir = sample.bt2_db_dir
            bt2_db_name = "pangenomes"

            if args.species_list:
                species_profile = sample.select_species(args.genome_coverage, args.species_list)
            else:
                species_profile = sample.select_species(args.genome_coverage)

            local_toc = download_reference(outputs.genomes, bt2_db_dir)
            db = UHGG(local_toc)

            # Because the pan genes are all named the same, so we NEED the species
            # subdirectories
            centroids_files = db.fetch_centroids(species_profile.keys(), bt2_db_dir)
            # Download centroids.ffn for every species in the restricted species profile.
            build_bowtie2_db(tempdir, bt2_db_name, centroids_files)
        # Perhaps avoid this giant conglomerated file, fetching instead submaps for each species.
        # TODO: Also colocate/cache/download in master for multiple slave subcommand invocations.

        bowtie2_align(args, bt2_db_dir, bt2_db_name, sort_aln=False)

        # Compute coverage of pangenome for each present species and write results to disk
        marker_genes_map = f"{marker_genes}/phyeco.map.lz4"
        species, genes = scan_centroids(centroids_files)
        num_covered_genes, species_mean_coverage, covered_genes = count_mapped_bp(args, tempdir, genes)
        markers = scan_markers(genes, marker_genes_map)
        species_markers_coverage = normalize(genes, covered_genes, markers)

        write_results(args.midas_outdir, species, num_covered_genes, species_markers_coverage, species_mean_coverage)
    except:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error.  Specify --debug flag to keep.")
            command(f"rm -rf {tempdir}", check=False)


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_run_genes(args)
