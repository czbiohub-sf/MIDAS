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
from iggtools.params.schemas import genes_summary_schema, genes_coverage_schema, MARKER_INFO_SCHEMA, format_data
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
    subparser.add_argument('--species_profile_path',
                           dest='species_profile_path',
                           type=str,
                           metavar="CHAR",
                           help=f"Species profile path for the prebuild bowtie2 index")

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


def scan_centroids(centroids_file, species_id):
    """ parse centroid fasta for given species_id"""
    centroids = {}
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


def write_genes_summary(genes_summary_list, genes_summary_path):
    with OutputStream(genes_summary_path) as stream:
        stream.write("\t".join(genes_summary_schema.keys()) + "\n")

        for record in genes_summary_list:
            stream.write("\t".join(map(format_data, record)) + "\n")


def write_genes_coverage(centroids, coverage_path):

    with OutputStream(coverage_path) as stream:
        stream.write('\t'.join(genes_coverage_schema.keys()) + '\n')

        for gene_id, gene_dict in centroids.items():
            if gene_dict["depth"] == 0: # Sparse by default here.
                continue
            vals = [gene_id, gene_dict["aligned_reads"], gene_dict["mapped_reads"], gene_dict["depth"], gene_dict["copies"]]
            stream.write("\t".join(map(format_data, vals)) + "\n")
    return "success"


def compute_gene_coverage(packed_args):
    """ Count number of bp mapped to each pan-gene. """
    pangenome_bamfile, gene_id, gene_length = packed_args

    with AlignmentFile(pangenome_bamfile) as bamfile:
        aligned_reads = bamfile.count(gene_id)
        mapped_reads = bamfile.count(gene_id, read_callback=keep_read)
        depth = sum((len(aln.query_alignment_sequence) / gene_length for aln in bamfile.fetch(gene_id)))
    return (aligned_reads, mapped_reads, depth)


def compute_median_marker_depth(centroids, species_id):
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


def compute_species_coverage(species_id, centroids_file, pangenome_bamfile, coverage_path):
    """ Compute coverage of pangenome for given species_id and write results to disk """
    global global_args
    args = global_args

    centroids = scan_centroids(centroids_file, species_id)
    #centroids = {k: centroids[k] for k in list(centroids)[:100]}

    args_list = []
    for gene_id in centroids.keys():
        args_list.append((pangenome_bamfile, gene_id, centroids[gene_id]["length"]))
    results = multiprocessing_map(compute_gene_coverage, args_list, num_procs=num_physical_cores)

    for gene_index, gene_id in enumerate(centroids):
        gene = centroids.get(gene_id)
        gene["aligned_reads"] = results[gene_index][0]
        gene["mapped_reads"] = results[gene_index][1]
        gene["depth"] = results[gene_index][2]

    marker_median_depth = compute_median_marker_depth(centroids, species_id)

    # Infer gene copy count
    if marker_median_depth > 0:
        for gene in centroids.values():
            gene["copies"] = gene["depth"] / marker_median_depth
    else:
        tsprint("Caution: the marker_gene_median_depth for {species_id} is zero")

    # Write gene coverage to file
    flag = write_genes_coverage(centroids, coverage_path)
    assert flag == "success", f"Error for write genes coverage for {species_id}"

    # Collect species-level stats
    nz_gene_depth = [gd["depth"] for gd in centroids.values() if gd["depth"] > 0]
    num_covered_genes = len(nz_gene_depth) #count # covered genes
    mean_coverage = np.mean(nz_gene_depth) # average gene depth
    pangenome_size = len(centroids)
    fraction_covered = num_covered_genes / pangenome_size
    aligned_reads = sum(g["aligned_reads"] for g in centroids.values())
    mapped_reads = sum(g["mapped_reads"] for g in centroids.values())

    values = [species_id, pangenome_size, num_covered_genes, fraction_covered, \
              mean_coverage, marker_median_depth, aligned_reads, mapped_reads]

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
            assert bowtie2_index_exists(bt2_db_dir, bt2_db_name) and os.path.exists(args.species_profile), f"Check the path bt2_db_dir and exists of species_profile_path"
            tsprint("Prebuild bowtie2 index and species_profile exit. Use them")

            bt2_db_dir = os.path.dirname(args.bt2_db_indexes)
            bt2_db_name = os.path.basename(args.bt2_db_indexes)

            species_profile = {}
            with InputStream() as stream:
                for species_id, coverage in select_from_tsv(stream, ["species_id", "coverage"]):
                    species_profile[species_id] = coverage

            sample.create_species_subdir(species_profile.keys(), args.debug, "dbs")
        else:
            bt2_db_dir = sample.get_target_layout("dbsdir")
            bt2_db_name = "pangenomes"
            bt2_db_temp_dir = sample.get_target_layout("dbs_tempdir")

            if args.species_list:
                species_profile = sample.select_species(args.genome_coverage, args.species_list)
            else:
                species_profile = sample.select_species(args.genome_coverage)

            sample.create_species_subdir(species_profile.keys(), args.debug, "dbs")

            local_toc = download_reference(outputs.genomes, bt2_db_dir)
            db = UHGG(local_toc)
            centroids_files = db.fetch_centroids(species_profile.keys(), bt2_db_temp_dir)
            build_bowtie2_db(bt2_db_dir, bt2_db_name, centroids_files)
            # Perhaps avoid this giant conglomerated file, fetching instead submaps for each species.
            # TODO: Also colocate/cache/download in master for multiple slave subcommand invocations.

        pangenome_bamfile = sample.get_target_layout("genes_pangenomes_bam")
        bowtie2_align(bt2_db_dir, bt2_db_name, pangenome_bamfile, args)
        samtools_index(pangenome_bamfile, args.debug)

        genes_summary_list = []
        for species_index, species_id in enumerate(list(species_profile.keys())):
            coverage_path = sample.get_target_layout("genes_coverage", species_id)
            genes_summary_list.append(compute_species_coverage(species_id, centroids_files[species_index], pangenome_bamfile, coverage_path))

        write_genes_summary(genes_summary_list, sample.get_target_layout("genes_summary"))

    except:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error.  Specify --debug flag to keep.")
            sample.remove_output_dir()
        raise


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_run_genes(args)
