import json
import os
import multiprocessing

from collections import defaultdict
import numpy as np
import Bio.SeqIO
from pysam import AlignmentFile  # pylint: disable=no-name-in-module

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, command, InputStream, OutputStream, select_from_tsv, multithreading_map, multiprocessing_map, download_reference, num_physical_cores
from iggtools.params import outputs
from iggtools.common.bowtie2 import build_bowtie2_db, bowtie2_align, samtools_index, bowtie2_index_exists
from iggtools.models.uhgg import UHGG, pangenome_file, marker_genes_mapfile
from iggtools.params.schemas import genes_profile_schema, genes_info_schema, genes_schema, MARKER_INFO_SCHEMA, format_data
from iggtools.models.sample import Sample


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
    subparser.add_argument('--sample_name',
                           dest='sample_name',
                           required=True,
                           help="Unique sample identifier")
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
    subparser.add_argument('--aln_sort',
                           action='store_true',
                           default=True,
                           help=f"Sort BAM file.")

    subparser.add_argument('--max_reads',
                           dest='max_reads',
                           type=int,
                           metavar="INT",
                           help=f"Number of reads to use from input file(s).  (All)")
    return main_func


def count_mapped_bp(pangenome_bamfile, genes):
    """ Count number of bp mapped to each gene across pangenomes.
    Return number covered genes and average gene depth per species.
    Result contains only covered species, but being a defaultdict,
    would yield 0 for any uncovered species, which is appropriate.
    """

    global global_args
    args = global_args

    bamfile = AlignmentFile(pangenome_bamfile, "rb")
    covered_genes = {}

    # loop over alignments, sum values per gene
    for aln in bamfile.fetch(until_eof=True):
        gene_id = bamfile.get_reference_name(aln.reference_id)
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
        if gene_depth > 0:  # WRONG. covered_genes is the same with centriods. This should always pass, because args.aln_cov is always >0.
            species_id = g["species_id"]
            nonzero_gene_depths[species_id].append(gene_depth)

    # Compute number of covered genes per species, and average gene depth.
    num_covered_genes = defaultdict(int)
    mean_coverage = defaultdict(float)
    for species_id, non_zero_depths in nonzero_gene_depths.items():
        num_covered_genes[species_id] = len(non_zero_depths)
        mean_coverage[species_id] = np.mean(non_zero_depths)

    return num_covered_genes, mean_coverage, covered_genes



def write_genes_coverage(my_args):
    species_genes, path = my_args
    """ Write centroid genes mappign for one species """
    with OutputStream(path) as sp_out:
        sp_out.write('\t'.join(genes_schema.keys()) + '\n')
        for gene_id, gene in species_genes.items():
            if gene["depth"] == 0:
                # Sparse by default here.
                # You can get the pangenome_size from the summary file, emitted below.
                continue
            values = [gene_id, gene["mapped_reads"], gene["depth"], gene["copies"]]
            sp_out.write("\t".join(map(format_data, values)) + "\n")


def write_genes_summary(species, num_covered_genes, species_markers_coverage, species_mean_coverage, outfile):
    # summary stats
    with OutputStream(outfile) as file:
        file.write('\t'.join(genes_profile_schema.keys()) + '\n')
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


def scan_centroids(centroids_file):
    centroids = {}
    # Read in the per-species centroid fasta file
    # TODO part of the information should already be processed during the database build
    with InputStream(centroids_file) as file:
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
            centroids[centroid_gene_id] = centroid_gene
    return centroids


def keep_read(aln):
    global global_args
    args = global_args

    align_len = len(aln.query_alignment_sequence)
    query_len = aln.query_length
    # min pid
    if 100 * (align_len - dict(aln.tags)['NM']) / float(align_len) < args.aln_mapid:
        return False
    # min read quality
    if np.mean(aln.query_qualities) < args.aln_readq:
        return False
    # min map quality
    if aln.mapping_quality < args.aln_mapq:
        return False
    # min aln cov
    if align_len / float(query_len) < args.aln_cov:
        return False
    return True


def compute_gene_coverage(packed_args):
    pangenome_bamfile, gene_id, gene_length = packed_args
    with AlignmentFile(pangenome_bamfile) as bamfile:
        aligned_reads = bamfile.count(gene_id)
        mapped_reads = bamfile.count(gene_id, read_callback=keep_read)
        depth = sum((len(aln.query_alignment_sequence) / gene_length for aln in bamfile.fetch(gene_id)))
        return (aligned_reads, mapped_reads, depth)


def compute_median_marker_depth(species_id, centroids):
    """ Compute median marker depth to infer gene copy count """

    # Get the gene_id - marker_id map
    awk_command = f"awk \'$1 == \"{species_id}\"\'"
    markers = dict()
    with InputStream(marker_genes_mapfile(), awk_command) as stream:
        for gene_id, marker_id in select_from_tsv(stream, ["gene_id", "marker_id"], schema=MARKER_INFO_SCHEMA):
            assert marker_id not in markers, f"Interesting, marker_id {marker_id} for species {species_id} corresponds to multiple gene_ids."
            markers[gene_id] = marker_id

    # Get the gene_id to centroid_gene_id map
    gene_info = dict()
    with InputStream(pangenome_file(species_id, "gene_info.txt.lz4")) as stream:
        for gene_id, centroids_gene_id in select_from_tsv(stream, ["gene_id", "centroid_99"]):
            if gene_id in markers.keys():
                gene_info[gene_id] = centroids_gene_id

    # Normalize gene_depth by median_marker_depth
    markers_depth = defaultdict(int)
    for gene_id in markers.keys():
        marker_id = markers[gene_id]
        centroid_gene_id = gene_info[gene_id]
        markers_depth[marker_id] += centroids[centroid_gene_id]["depth"]

    # Compute median marker depth for current species_id
    return np.median(list(markers_depth.values()))


def compute_species_coverage(species_id, centroids_file, pangenome_bamfile, path):

    global global_args
    args = global_args

    centroids = scan_centroids(centroids_file)
    centroids = {k: centroids[k] for k in list(centroids)[:500]}

    args_list = []
    for gene_id in centroids.keys():
        args_list.append((pangenome_bamfile, gene_id, centroids[gene_id]["length"]))
    results = multiprocessing_map(compute_gene_coverage, args_list, num_procs=num_physical_cores)

    for gene_index, gene_id in enumerate(centroids):
        gene = centroids.get(gene_id)
        gene["aligned_reads"] = results[gene_index][0]
        gene["mapped_reads"] = results[gene_index][1]
        gene["depth"] = results[gene_index][2]

    # Calculate the genes covered with non-zero depth
    nz_gene_depth = [gd["depth"] for gd in centroids.values() if gd["depth"] > 0]
    num_covered_genes = len(nz_gene_depth)
    mean_coverage = np.mean(nz_gene_depth)
    print(mean_coverage, num_covered_genes)

    marker_median_depth = compute_median_marker_depth(species_id, centroids)
    # Infer gene copy count
    if marker_median_depth > 0:
        for gene in centroids.values():
            gene["copies"] = gene["depth"] / marker_median_depth
    else:
        tsprint("Caution: the marker_gene_median_depth for {species_id} is zero")

    # write to file
    with OutputStream(path) as sp_out:
        sp_out.write('\t'.join(genes_schema.keys()) + '\n')
        for gene_id, gene in covered_genes.items():
            values = [gene_id, gene["mapped_reads"], gene["depth"], gene["copies"]]
            sp_out.write("\t".join(map(format_data, values)) + "\n")

    # summary
    pangenome_size = len(centroids)
    aligned_reads = sum(g["aligned_reads"] for g in covered_genes.values())
    mapped_reads = sum(g["mapped_reads"] for g in covered_genes.values())
    values = [species_id, pangenome_size, \
              num_covered_genes, num_covered_genes / pangenome_size, \
              mean_coverage, markers_coverage, aligned_reads, mapped_reads]
    print(values)
    return values


def midas_run_genes(args):

    sample = Sample(args.sample_name, args.midas_outdir, "genes")
    sample.create_output_dir(args.debug)

    global global_sample
    global_sample = sample
    global global_args
    global_args = args

    try:
        if args.bt2_db_indexes:
            if bowtie2_index_exists(os.path.dirname(args.bt2_db_indexes), os.path.basename(args.bt2_db_indexes)):
                print("good the pre built bt2 index exist")
                bt2_db_dir = os.path.dirname(args.bt2_db_indexes)
                bt2_db_name = os.path.basename(args.bt2_db_indexes)
            else:
                print("Error: good the pangenome pre built bt2 index exist")
                exit(0)
        else:
            bt2_db_dir = sample.get_target_layout("dbsdir")
            bt2_db_temp_dir = sample.get_target_layout("dbs_tempdir")
            bt2_db_name = "pangenomes"

            if args.species_list:
                species_profile = sample.select_species(args.genome_coverage, args.species_list)
            else:
                species_profile = sample.select_species(args.genome_coverage)

            local_toc = download_reference(outputs.genomes, bt2_db_dir)
            db = UHGG(local_toc)

            sample.create_species_subdir(species_profile.keys())
            centroids_files = db.fetch_centroids(species_profile.keys(), bt2_db_temp_dir)

            build_bowtie2_db(bt2_db_dir, bt2_db_name, centroids_files)
            # Perhaps avoid this giant conglomerated file, fetching instead submaps for each species.
            # TODO: Also colocate/cache/download in master for multiple slave subcommand invocations.

        # Here: database and bam file should not have the default behavior as same path ...
        # This doesn't make sense.
        bowtie2_align(bt2_db_dir, bt2_db_name, args)
        samtools_index(args, bt2_db_dir, bt2_db_name)

        # TODO: the following compute is not paralled. Finish the workflow for one species,
        # compute and write for one species. Then do multiprocessing.

        if False:
            # Compute coverage of pangenome for each present species and write results to disk
            marker_genes_map = marker_genes_mapfile()

            print("scan_centroids")
            species, genes = scan_centroids(species_profile.keys(), centroids_files)
            print("count_mapped_bp")
            num_covered_genes, species_mean_coverage, covered_genes = count_mapped_bp(sample.get_target_layout("genes_pangenomes_bam"), genes)
            markers = scan_markers(genes, marker_genes_map)
            print("normalize")
            species_markers_coverage = normalize(genes, covered_genes, markers)
            # Write results
            print("write results")
            argument_list = []
            for species_id, species_genes in species.items():
                path = sample.get_target_layout("genes_coverage", species_id)
                argument_list.append((species_genes, path))
            contigs_files = multithreading_map(write_genes_coverage, argument_list, num_physical_cores)

            write_genes_summary(species, num_covered_genes, species_markers_coverage, species_mean_coverage, sample.get_target_layout("genes_summary"))

        species_ids = list(species_profile.keys())[:1]
        for species_index, species_id in enumerate(species_ids):
            pangenome_bamfile = sample.get_target_layout("genes_pangenomes_bam")
            coverage_path = sample.get_target_layout("genes_coverage", species_id)
            compute_species_coverage(species_id, centroids_files[species_index], pangenome_bamfile, coverage_path)
    except:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error.  Specify --debug flag to keep.")
            sample.remove_output_dir()


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_run_genes(args)
