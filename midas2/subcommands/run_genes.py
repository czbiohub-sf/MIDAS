#!/usr/bin/env python3
import json
import os

from collections import defaultdict
from itertools import repeat
import numpy as np
from pysam import AlignmentFile  # pylint: disable=no-name-in-module

from midas2.common.argparser import add_subcommand
from midas2.common.utils import tsprint, InputStream, OutputStream, select_from_tsv, multiprocessing_map, args_string
from midas2.common.utilities import scan_cluster_info, scan_centroid_prev
from midas2.models.midasdb import MIDAS_DB
from midas2.models.sample import Sample
from midas2.models.species import Species, parse_species
from midas2.params.schemas import genes_summary_schema, genes_coverage_schema, format_data, DECIMALS6
from midas2.common.bowtie2 import build_bowtie2_db, bowtie2_align, samtools_index, bowtie2_index_exists, _keep_read
from midas2.params.inputs import MIDASDB_NAMES


DEFAULT_ALN_MAPID = 94.0
DEFAULT_ALN_MAPQ = 0
DEFAULT_ALN_COV = 0.75
DEFAULT_ALN_READQ = 20

DEFAULT_READ_DEPTH = 0
DEFAULT_MARKER_MEDIAN_DEPTH = 2.0

DEFAULT_MAX_FRAGLEN = 500
DEFAULT_NUM_CORES = 8


def register_args(main_func):
    subparser = add_subcommand('run_genes', main_func, help='Metagenomic pan-genome profiling')

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

    # Prebuilt/predownload
    subparser.add_argument('--prebuilt_bowtie2_indexes',
                           dest='prebuilt_bowtie2_indexes',
                           type=str,
                           metavar="CHAR",
                           help="Path to prebuilt pan-genome bowtie2 database indexes")
    subparser.add_argument('--prebuilt_bowtie2_species',
                           dest='prebuilt_bowtie2_species',
                           type=str,
                           metavar="CHAR",
                           help="Path to list of species in the prebuild bowtie2 indexes.")
    subparser.add_argument('--midasdb_name',
                           dest='midasdb_name',
                           type=str,
                           default="uhgg",
                           choices=MIDASDB_NAMES,
                           help="MIDAS Database name.")
    subparser.add_argument('--midasdb_dir',
                           dest='midasdb_dir',
                           type=str,
                           default="midasdb",
                           help="Path to local MIDAS Database.")

    subparser.add_argument('--species_list',
                           dest='species_list',
                           type=str,
                           metavar="CHAR",
                           help="Comma separated list of species ids OR path to list of species TXT.")
    subparser.add_argument('--select_by',
                           dest='select_by',
                           type=str,
                           default="median_marker_coverage",
                           help="Comma separated columns from species_profile to filter species.")
    subparser.add_argument('--select_threshold',
                           dest='select_threshold',
                           type=str,
                           metavar="CHAR",
                           default=str(DEFAULT_MARKER_MEDIAN_DEPTH),
                           help=f"Comman separated correponsding cutoff to select_by (>XX) ({DEFAULT_MARKER_MEDIAN_DEPTH}, )")

    #  Alignment flags (Bowtie2, or postprocessing)
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
                           help='Global/local read alignment (Default local, corresponds to the bowtie2 --local; global corresponds to the bowtie2 --end-to-end).')
    subparser.add_argument('--aln_interleaved',
                           action='store_true',
                           default=False,
                           help='FASTA/FASTQ file in -1 are paired and contain forward AND reverse reads')
    subparser.add_argument('--fragment_length',
                           type=float,
                           dest='fragment_length',
                           metavar='FLOAT',
                           default=DEFAULT_MAX_FRAGLEN,
                           help="Maximum fragment length for paired reads.")
    subparser.add_argument('--max_reads',
                           dest='max_reads',
                           type=int,
                           metavar="INT",
                           help="Number of reads to use from input file(s) for read alignment.  (All)")
    subparser.add_argument('--aln_extra_flags',
                           type=str,
                           dest='aln_extra_flags',
                           default='',
                           help='Extra bowtei2 align flags. E.g. --mm --ignore-quals')


    # Coverage flags
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
    subparser.add_argument('--aln_readq',
                           dest='aln_readq',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_ALN_READQ,
                           help=f"Discard reads with mean quality < READQ ({DEFAULT_ALN_READQ})")
    subparser.add_argument('--aln_cov',
                           dest='aln_cov',
                           default=DEFAULT_ALN_COV,
                           type=float,
                           metavar="FLOAT",
                           help=f"Discard reads with alignment coverage < ALN_COV ({DEFAULT_ALN_COV}).  Values between 0-1 accepted.")

    # Gene filters
    subparser.add_argument('--read_depth',
                           dest='read_depth',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_READ_DEPTH,
                           help=f"Discard genes with post-filtered reads < READ_DEPTH  ({DEFAULT_READ_DEPTH})")

    # Resource related
    subparser.add_argument('--num_cores',
                           dest='num_cores',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_NUM_CORES,
                           help=f"Number of physical cores to use ({DEFAULT_NUM_CORES})")
    # Prune centroids_99 for cleaner reads mapping
    subparser.add_argument('--prune_centroids99',
                           action='store_true',
                           default=False,
                           help='Prune shorter centroids99 within each centroids95 cluster')
    subparser.add_argument('--prune_opts',
                           dest='prune_opts',
                           type=str,
                           default="0.4",
                           help="Prune centroids99 shorter than 40% of its centroids.95.")
    subparser.add_argument('--aln_only',
                           dest = 'aln_only',
                           action='store_true',
                           default=False,
                           help='Perform only alignment steps.')
    subparser.add_argument('--keep_bam',
                           dest = 'keep_bam',
                           action='store_true',
                           default=False,
                           help='Keep BAM file.')
    return main_func


def keep_read(aln):
    global global_args
    args = global_args
    return _keep_read(aln, args.aln_mapid, args.aln_readq, args.aln_mapq, args.aln_cov)


def prepare_species(species_to_analyze, midas_db):
    global dict_of_species
    dict_of_species = {species_id: Species(species_id) for species_id in species_to_analyze}

    for species_id in species_to_analyze:
        sp = dict_of_species[species_id]
        sp.get_cluster_info_fp(midas_db)
        sp.get_c95_prevalence_fp(midas_db)

    return True


def collect_reads(c99_info, c95_prevalence, pangenome_bamfile):
    list_of_c99_ids = sorted(list(c99_info.keys()))

    c95_values = defaultdict(dict)
    with AlignmentFile(pangenome_bamfile) as bamfile:
        # Competitive alignment is done on centroid_99 level, while coverage are computed on centroid_95 level.
        for c99_id in list_of_c99_ids:
            c95_id = c99_info[c99_id]["centroid_95"]

            # centroids_95 is subset of centroids_99
            c99_length = c99_info[c99_id]["centroid_99_length"]
            c95_length = c99_info[c95_id]["centroid_99_length"]

            aligned_reads = bamfile.count(c99_id)
            mapped_reads = bamfile.count(c99_id, read_callback=keep_read)

            # Compute total per-position depth for aligned gene region
            c99_covered_bases = 0
            c99_total_depth = 0
            counts = bamfile.count_coverage(c99_id, read_callback=keep_read)
            for within_chunk_index in range(0, c99_length):
                # Per-position depth: total number of bases mappped to given c99
                gene_depth = sum([counts[nt][within_chunk_index] for nt in range(4)])
                c99_total_depth += gene_depth
                if gene_depth > 0:
                    c99_covered_bases += 1

            if c99_total_depth == 0: # Sparse by default.
                continue

            if c95_id not in c95_values:
                c95_values[c95_id]["c95_id"] = c95_id
                c95_values[c95_id]["c95_length"] = c95_length
                c95_values[c95_id]["aligned_reads"] = aligned_reads
                c95_values[c95_id]["mapped_reads"] = mapped_reads
                c95_values[c95_id]["read_depth"] = c99_total_depth #read depths are summed at c95
                c95_values[c95_id]["mean_coverage"] = 0.0
                c95_values[c95_id]["copy_number"] = 0.0
                c95_values[c95_id]["c95_prevalence"] = c95_prevalence[c95_id]
                c95_values[c95_id]["marker_id"] = c99_info[c95_id]["centroid_95_marker_id"] if c99_info[c95_id]["is_centroid_95_marker"] else ""
            else:
                c95_values[c95_id]["aligned_reads"] += aligned_reads
                c95_values[c95_id]["mapped_reads"] += mapped_reads
                c95_values[c95_id]["read_depth"] += c99_total_depth
    return c95_values


def decorate_coverage(c95_values):
    # Compute the median read coverage for all mapped marker genes for given species
    list_of_markers = set(r['marker_id'] for r in c95_values.values() if r['marker_id'])
    markers_coverage = dict(zip(list_of_markers, [0.0]*len(list_of_markers)))
    for c95_cov in c95_values.values():
        if c95_cov["is_centroid_95_marker"]:
            markers_coverage[c95_cov["marker_id"]] += c95_cov["read_depth"] / c95_cov["c95_length"]
    median_marker_coverage = np.median(list(map(lambda x: float(format(x, DECIMALS6)), markers_coverage.values())))

    # Update the mean_coverage and copy_number values of c95_values
    pangenome_size = 0
    for c95_cov in c95_values.values():
        c95_cov["mean_coverage"] = float(c95_cov["read_depth"] / c95_cov["gene_length"])
        if median_marker_coverage > 0:
            c95_cov["copy_number"] = float(c95_cov["mean_coverage"] / median_marker_coverage)
        pangenome_size += 1

    return median_marker_coverage, pangenome_size


def process_species(species_id):
    """ Collect total number of read depths over each covered position for centroids_99 """

    global global_args
    global sample
    pangenome_bamfile = sample.get_target_layout("genes_pangenomes_bam")

    global dict_of_species
    sp = dict_of_species[species_id]

    c99_info = scan_cluster_info(sp.cluster_info_fp, "centroid_99")
    c95_prevalence = scan_centroid_prev(sp.c95_prevalence_fp, "95")

    c95_values = collect_reads(c99_info, c95_prevalence, pangenome_bamfile)
    median_marker_coverage, pangenome_size = decorate_coverage(c95_values)

    # Aggreate species coverage summary AND write per-gene coverage to file
    c95_summary = {
        "species_id": species_id,
        "pangenoem_size": pangenome_size,
        "covered_genes": 0,
        "fraction_covered": 0,
        "aligned_reads": 0,
        "mapped_reads": 0,
        "mean_coverage": 0,
        "marker_coverage": median_marker_coverage,
    }

    c95_cov_fp = sample.get_target_layout("genes_coverage", species_id)
    with OutputStream(c95_cov_fp) as stream:
        stream.write('\t'.join(genes_coverage_schema.keys()) + '\n')
        for rec in c95_values.values():
            if rec["mapped_reads"] > global_args.read_depth:
                # accumuate across covered centroid_95s
                c95_summary["covered_genes"] += 1
                c95_summary["aligned_reads"] += rec["aligned_reads"]
                c95_summary["mapped_reads"] += rec["mapped_reads"]
                c95_summary["mean_coverage"] += rec["mean_coverage"]
                # write single centroid_95 coverage to file
                stream.write("\t".join(map(format_data, rec, repeat(DECIMALS6, len(rec)))) + "\n")

    c95_summary["fraction_covered"] = c95_summary["covered_genes"] / c95_summary["pangenome_size"]
    c95_summary["mean_coverage"] = c95_summary["mean_coverage"] / c95_summary["covered_genes"] if c95_summary["covered_genes"] > 0 else 0
    return c95_summary


def write_species_summary(species_covs, genes_coverage_path):
    with OutputStream(genes_coverage_path) as stream:
        stream.write("\t".join(genes_summary_schema.keys()) + "\n")
        for rec in species_covs:
            stream.write("\t".join(map(format_data, rec.values())) + "\n")


def run_genes(args):

    try:
        global global_args
        global_args = args

        global sample
        sample = Sample(args.sample_name, args.midas_outdir, "genes")
        sample.create_dirs(["outdir", "tempdir"], args.debug)

        with OutputStream(sample.get_target_layout("genes_log")) as stream:
            stream.write(f"Single sample pan-gene copy number variant calling in subcommand {args.subcommand} with args\n{json.dumps(args_string(args), indent=4)}\n")

        species_list = parse_species(args)

        if args.prebuilt_bowtie2_indexes:
            bt2_db_dir = os.path.dirname(args.prebuilt_bowtie2_indexes)
            bt2_db_name = os.path.basename(args.prebuilt_bowtie2_indexes)

            assert bowtie2_index_exists(bt2_db_dir, bt2_db_name), f"Provided {bt2_db_dir}/{bt2_db_name} don't exist."
            assert (args.prebuilt_bowtie2_species and os.path.exists(args.prebuilt_bowtie2_species)), "Require list of speices used to build the provided Bowtie2 indexes."

            tsprint(f"Read in list of species used to build provided bowtie2 indexes {bt2_db_dir}/{bt2_db_name}")
            with InputStream(args.prebuilt_bowtie2_species) as stream:
                bt2_species_list = [spid[0] for spid in select_from_tsv(stream, schema={"species_id": str})]

            # Update species_list: either/or
            species_list = list(set(species_list) & set(bt2_species_list)) if species_list else bt2_species_list
        else:
            sample.create_dirs(["bt2_indexes_dir"], args.debug)
            bt2_db_dir = sample.get_target_layout("bt2_indexes_dir")
            bt2_db_name = "pangenomes"

        # Select abundant species present in the sample for SNPs calling
        select_thresholds = args.select_threshold.split(',')
        no_filter = len(select_thresholds) == 1 and float(select_thresholds[0]) == -1

        species_to_analyze = species_list if no_filter else sample.select_species(args, species_list)
        species_counts = len(species_to_analyze)

        assert species_counts > 0, "No (specified) species pass the marker_depth filter, please adjust the marker_depth or species_list"
        tsprint(len(species_to_analyze))

        tsprint("MIDAS2::prepare_species::start")
        num_cores = min(species_counts, args.num_cores)
        midas_db = MIDAS_DB(os.path.abspath(args.midasdb_dir), args.midasdb_name, num_cores)
        midas_db.fetch_files("pangenome", species_to_analyze)
        prepare_species(species_to_analyze, midas_db, args.chunk_size)
        tsprint("MIDAS2::prepare_species::finish")

        # Build Bowtie indexes for species in the restricted species profile
        tsprint("MIDAS2::build_bowtie2db::start")
        if args.prune_centroids:
            centroids_files = midas_db.fetch_files("pruned_centroids", species_to_analyze)
        else:
            centroids_files = midas_db.fetch_files("pangenome_centroids", species_to_analyze)
        build_bowtie2_db(bt2_db_dir, bt2_db_name, centroids_files, args.num_cores)
        tsprint("MIDAS2::build_bowtie2db::finish")

        # Align reads to pangenome database
        tsprint("MIDAS2::bowtie2_align::start")
        sample.create_species_subdirs(species_to_analyze, "temp", args.debug)
        pangenome_bamfile = sample.get_target_layout("genes_pangenomes_bam")
        bowtie2_align(bt2_db_dir, bt2_db_name, pangenome_bamfile, args)
        samtools_index(pangenome_bamfile, args.debug, args.num_cores)
        tsprint("MIDAS2::bowtie2_align::finish")

        if args.aln_only:
            return

        tsprint("MIDAS2::multiprocessing_map::start")
        species_coverages = multiprocessing_map(process_species, species_to_analyze, num_cores)
        tsprint("MIDAS2::multiprocessing_map::finish")

        tsprint("MIDAS2::write_species_summary::start")
        write_species_summary(species_coverages, sample.get_target_layout("genes_summary"))
        tsprint("MIDAS2::write_species_summary::finish")

        if not args.keep_bam:
            sample.remove_dirs(["tempdir"])
        if not args.prebuilt_bowtie2_indexes:
            sample.remove_dirs(["bt2_indexes_dir"])

    except Exception as error:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error.  Specify --debug flag to keep.")
            sample.remove_dirs(["outdir", "tempdir"])
            if not args.prebuilt_bowtie2_indexes:
                sample.remove_dirs(["bt2_indexes_dir"])
        raise error


@register_args
def main(args):
    tsprint(f"Single sample pan-gene copy number variant calling in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    run_genes(args)
