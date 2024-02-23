#!/usr/bin/env python3
import json
import os
from math import ceil
from collections import defaultdict
from itertools import repeat
import numpy as np
from pysam import AlignmentFile  # pylint: disable=no-name-in-module

from midas2.common.argparser import add_subcommand
from midas2.common.utils import tsprint, InputStream, OutputStream, select_from_tsv, multiprocessing_map, args_string, command, multithreading_map
from midas2.common.utilities import extract_genomeid
from midas2.models.midasdb import MIDAS_DB
from midas2.models.sample import Sample
from midas2.models.species import Species, parse_species
from midas2.params.schemas import genes_summary_schema, fetch_genes_depth_schema, format_data, DECIMALS6, fetch_genes_chunk_schema
from midas2.common.bowtie2 import build_bowtie2_db, bowtie2_align, samtools_index, samtools_idxstats, bowtie2_index_exists, _keep_read
from midas2.params.inputs import MIDASDB_NAMES


DEFAULT_ALN_MAPID = 94.0
DEFAULT_ALN_MAPQ = 0
DEFAULT_ALN_COV = 0.75
DEFAULT_ALN_READQ = 20

DEFAULT_TOTAL_DEPTH = 1
DEFAULT_MARKER_MEDIAN_DEPTH = 2.0

DEFAULT_MAX_FRAGLEN = 500
DEFAULT_NUM_CORES = 8
DEFAULT_PRUNE_CUTOFF = 0.4
DEFAULT_CLUSTER_ID = '99'
DEFAULT_MARKER_CUTOFF = 0.01


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
                           default="median_marker_depth",
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
                           default='global',
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
                           metavar="CHAR",
                           dest='aln_extra_flags',
                           default='',
                           help="Extra bowtei2 align flags. E.g. --mm --ignore-quals")


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
    subparser.add_argument('--total_depth',
                           dest='total_depth',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_TOTAL_DEPTH,
                           help=f"Discard genes with post-filtered reads < total_depth  ({DEFAULT_TOTAL_DEPTH})")

    # Resource related
    subparser.add_argument('--num_cores',
                           dest='num_cores',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_NUM_CORES,
                           help=f"Number of physical cores to use ({DEFAULT_NUM_CORES})")
    # Prune centroids_99 for cleaner reads mapping
    subparser.add_argument('--prune_centroids',
                           action='store_true',
                           default=False,
                           help='Prune shorter centroids99 within each centroids95 cluster')
    subparser.add_argument('--prune_method',
                           dest='prune_method',
                           type=str,
                           default='max',
                           choices=['max', 'median', 'mean'],
                           help="Prune methods: max, median or mean.")
    subparser.add_argument('--prune_cutoff',
                           dest='prune_cutoff',
                           type=float,
                           default=DEFAULT_PRUNE_CUTOFF,
                           help=f"Prune cutoff: for each centroid_95, centroid99 shorter than {DEFAULT_PRUNE_CUTOFF} of the chosen method are pruned for reading mapping.")
    subparser.add_argument('--remove_singleton',
                           action='store_true',
                           default=False,
                           help='Remove c75 clusters with only one gene member in species with more genomes.')

    subparser.add_argument('--alignment_only',
                           dest = 'alignment_only',
                           action='store_true',
                           default=False,
                           help='Perform only alignment steps.')
    subparser.add_argument('--cluster_level',
                           dest='cluster_level',
                           type=str,
                           default=DEFAULT_CLUSTER_ID,
                           choices=['75', '80', '85', '90', '95', '99'],
                           help=f"Aggregate reads to cluster level other than {DEFAULT_CLUSTER_ID}. Only for debug purpose.")

    subparser.add_argument('--remove_bam',
                           dest = 'remove_bam',
                           action='store_true',
                           default=False,
                           help='Remove BAM file.')
    subparser.add_argument('--remove_bt2_index',
                           dest = 'remove_bt2_index',
                           action='store_true',
                           default=False,
                           help='Remove bowtie2 index files.')
    return main_func


def keep_read(aln):
    global global_args
    args = global_args
    return _keep_read(aln, args.aln_mapid, args.aln_readq, args.aln_mapq, args.aln_cov)


def fetch_prebuilt_bt2db(args, species_list):
    bt2_db_dir = os.path.dirname(args.prebuilt_bowtie2_indexes)
    bt2_db_name = os.path.basename(args.prebuilt_bowtie2_indexes)

    assert bowtie2_index_exists(bt2_db_dir, bt2_db_name), f"Provided {bt2_db_dir}/{bt2_db_name} don't exist."
    assert (args.prebuilt_bowtie2_species and os.path.exists(args.prebuilt_bowtie2_species)), "Require list of speices used to build the provided Bowtie2 indexes."

    tsprint(f"Read in list of species used to build provided bowtie2 indexes {bt2_db_dir}/{bt2_db_name}")
    with InputStream(args.prebuilt_bowtie2_species) as stream:
        bt2_species_list = [spid[0] for spid in select_from_tsv(stream, schema={"species_id": str})]
    # Update species_list: either/or
    species_list = list(set(species_list) & set(bt2_species_list)) if species_list else bt2_species_list
    return bt2_db_dir, bt2_db_name, species_list


def get_bowtie2_indexes(args, species_list):
    if args.prebuilt_bowtie2_indexes:
        bt2_db_dir, bt2_db_name, species_list = fetch_prebuilt_bt2db(args, species_list)
    else:
        sample.create_dirs(["bt2_indexes_dir"], args.debug)
        bt2_db_dir = sample.get_target_layout("bt2_indexes_dir")
        bt2_db_name = "pangenomes"
    return bt2_db_dir, bt2_db_name, species_list


def filter_species_list(args, species_list):
    select_thresholds = args.select_threshold.split(',')
    no_filter = len(select_thresholds) == 1 and float(select_thresholds[0]) == -1

    species_to_analyze = species_list if no_filter else sample.select_species(args, species_list)
    return species_to_analyze


def fetch_pruned_centroids(midas_db, species_id, prune_method, prune_cutoff, remove_singleton):
    if remove_singleton:
        return midas_db.get_target_layout("pruned_centroids_rs", False, species_id, prune_method, prune_cutoff)
    return midas_db.get_target_layout("pruned_centroids", False, species_id, prune_method, prune_cutoff)


def _fetch_cxx_info(pargs):
    sp, midas_db, xx = pargs
    for xx in set(["99", xx]):
        sp.set_clusters_info_fp(midas_db, xx)
        sp.get_clusters_info(xx)
    return True


def prepare_species(species_to_analyze, midas_db):
    global dict_of_species
    global global_args
    dict_of_species = {species_id: Species(species_id) for species_id in species_to_analyze}
    num_cores = min(midas_db.num_cores, 8)
    multithreading_map(_fetch_cxx_info, [(sp, midas_db, global_args.cluster_level) for sp in dict_of_species.values()], num_cores)
    return True


def fetch_genes_from_bam(idxstats_fp, midas_db):
    """ Compute the genes present in the BAM file """
    species_for_genome = midas_db.uhgg.genomes
    global readonly_bamgenes
    readonly_bamgenes = {}
    with InputStream(idxstats_fp) as stream:
        for _ in select_from_tsv(stream, schema={"gene_id": str}):
            geneid = _[0]
            speciesid = species_for_genome[extract_genomeid(geneid)]
            readonly_bamgenes[geneid] = speciesid
    return True


def compute_pileup_per_chunk(pargs):
    """ Collect total number of read depths over each covered position for centroids_99 """
    global global_args
    global dict_of_species
    global readonly_bamgenes

    chunk_id, chunk_start, chunk_end, pangenome_bamfile, headerless_sliced_path, xx = pargs

    if global_args.debug and os.path.exists(headerless_sliced_path):
        tsprint(f"Skipping compute pileup for chunk {chunk_id} in debug mode as temporary data exists: {headerless_sliced_path}")
        return headerless_sliced_path

    chunk_geneids_list = list(readonly_bamgenes.keys())[chunk_start:chunk_end]

    cxx_values = defaultdict(dict)
    with AlignmentFile(pangenome_bamfile) as bamfile:
        # Competitive alignment is done on centroid_99 level.
        for c99_id in chunk_geneids_list:
            c99_aligned_reads = bamfile.count(c99_id)
            if c99_aligned_reads < global_args.total_depth:
                continue

            c99_mapped_reads = bamfile.count(c99_id, read_callback=keep_read)
            if c99_mapped_reads < global_args.total_depth:
                continue

            species_id = readonly_bamgenes[c99_id]
            sp = dict_of_species[species_id]

            c99_info = sp.clusters_info['99']
            cxx_info = sp.clusters_info[xx]

            cxx_id = c99_info[c99_id][f"centroid_{xx}"]
            c99_length = c99_info[c99_id]["centroid_99_gene_length"]
            cxx_length = cxx_info[cxx_id][f"centroid_{xx}_gene_length"]

            # Compute total per-position depth for aligned gene region
            c99_covered_bases = 0
            c99_total_depth = 0
            counts = bamfile.count_coverage(c99_id, read_callback=keep_read)
            for within_chunk_index in range(0, c99_length):
                # Per-position depth: total number of bases mappped
                gene_depth = sum([counts[nt][within_chunk_index] for nt in range(4)])
                c99_total_depth += gene_depth
                if gene_depth > 0:
                    c99_covered_bases += 1
            if c99_total_depth == 0: # Sparse by default.
                continue
            c99_mean_depth = float(c99_total_depth / c99_length)

            if cxx_id not in cxx_values:
                cxx_values[cxx_id] = {
                    "species_id": species_id,
                    f"c{xx}_id": cxx_id,
                    f"c{xx}_length": cxx_length,
                    "aligned_reads": c99_aligned_reads,
                    "mapped_reads": c99_mapped_reads,
                    "total_depth": c99_total_depth,
                    "mean_depth": c99_mean_depth,
                    "copy_number": 0.0,
                    "genome_prevalence": cxx_info[cxx_id][f"centroid_{xx}_genome_prevalence"],
                    "marker_id": cxx_info[cxx_id][f"centroid_{xx}_marker_id"],
                }
            else:
                cxx_values[cxx_id]["aligned_reads"] += c99_aligned_reads
                cxx_values[cxx_id]["mapped_reads"] += c99_mapped_reads
                cxx_values[cxx_id]["total_depth"] += c99_total_depth
                cxx_values[cxx_id]["mean_depth"] += c99_mean_depth

    with OutputStream(headerless_sliced_path) as stream:
        for rec in cxx_values.values():
            stream.write("\t".join(map(format_data, rec.values(), repeat(DECIMALS6, len(rec)))) + "\n")

    return headerless_sliced_path


def merge_depth_across_chunks(list_of_chunks_pileup, xx):
    depth_schema = fetch_genes_chunk_schema(xx)
    depth_cols = list(depth_schema.keys())[1:]

    # operational gene cluster level other than 99 could end up in different batches of c99s
    cxx_values = defaultdict(lambda: defaultdict(dict))
    for chunk_file in list_of_chunks_pileup:
        with InputStream(chunk_file) as stream:
            for r in select_from_tsv(stream, schema=depth_schema):
                species_id, cxx_id = r[:2]
                cxx_val = dict(zip(depth_cols, r[1:]))
                if cxx_id not in cxx_values[species_id]:
                    cxx_values[species_id][cxx_id] = cxx_val
                else:
                    cxx_values[species_id][cxx_id]["aligned_reads"] += cxx_val["aligned_reads"]
                    cxx_values[species_id][cxx_id]["mapped_reads"] += cxx_val["mapped_reads"]
                    cxx_values[species_id][cxx_id]["total_depth"] += cxx_val["total_depth"]
                    cxx_values[species_id][cxx_id]["mean_depth"] += cxx_val["mean_depth"]
    return cxx_values


def compute_median_marker_depth(cxx_values, list_of_markers):
    """ Compute the median read coverage for ALL marker genes for given species """
    markers_depth = dict(zip(list_of_markers, [0.0]*len(list_of_markers)))

    for cxx_cov in cxx_values.values():
        if cxx_cov["marker_id"]:
            markers_depth[cxx_cov["marker_id"]] += cxx_cov['mean_depth']

    median_marker_depth = np.median(list(map(lambda x: float(format(x, DECIMALS6)), markers_depth.values())))
    return median_marker_depth


def compute_species_summary(dict_of_gene_depth):
    global global_args
    global sample
    global dict_of_species

    xx = global_args.cluster_level
    species_depth_summary_dict = defaultdict(dict)
    for species_id, sp in dict_of_species.items():
        sp = dict_of_species[species_id]
        cxx_values = dict_of_gene_depth[species_id]

        cxx_info = sp.clusters_info[xx]
        list_of_markers = sp.list_of_markers[xx]
        median_marker_depth = compute_median_marker_depth(cxx_values, list_of_markers)
        tsprint(f"median marker depth for {species_id} is {median_marker_depth}")

        cxx_summary = {
            "species_id": species_id,
            "pangenome_size": len(cxx_info),
            "covered_genes": 0,
            "fraction_covered": 0,
            "aligned_reads": 0,
            "mapped_reads": 0,
            "mean_depth": 0,
            "marker_depth": median_marker_depth,
        }

        cxx_cov_fp = sample.get_target_layout("genes_depth", species_id)
        with OutputStream(cxx_cov_fp) as stream:
            stream.write('\t'.join(fetch_genes_depth_schema(xx).keys()) + '\n')
            for rec in cxx_values.values():
                # Collect species depth summary
                cxx_summary["covered_genes"] += 1
                cxx_summary["aligned_reads"] += rec["aligned_reads"]
                cxx_summary["mapped_reads"] += rec["mapped_reads"]
                cxx_summary["mean_depth"] += rec["mean_depth"]
                # Update copy_number in one pass
                rec["copy_number"] = float(rec["mean_depth"] / median_marker_depth) if median_marker_depth > 0 else 0.0
                # Write species, gene depth to file
                stream.write("\t".join(map(format_data, rec.values(), repeat(DECIMALS6, len(rec)))) + "\n")

        cxx_summary["fraction_covered"] = cxx_summary["covered_genes"] / cxx_summary["pangenome_size"]
        cxx_summary["mean_depth"] = cxx_summary["mean_depth"] / cxx_summary["covered_genes"] if cxx_summary["covered_genes"] > 0 else 0

        species_depth_summary_dict[species_id] = cxx_summary

    return species_depth_summary_dict


def write_species_summary(species_covs, genes_depth_path):
    with OutputStream(genes_depth_path) as stream:
        stream.write("\t".join(genes_summary_schema.keys()) + "\n")
        for rec in species_covs.values():
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

        # Read in list of species
        species_list = parse_species(args)
        bt2_db_dir, bt2_db_name, species_list = get_bowtie2_indexes(args, species_list)

        # Restricted species profile: only abundant species based on species SGC profiling results
        species_to_analyze = filter_species_list(args, species_list)
        species_counts = len(species_to_analyze)
        assert species_counts > 0, "No (specified) species pass the marker_depth filter, please adjust the marker_depth or species_list"
        tsprint(len(species_to_analyze))

        num_cores = min(species_counts, args.num_cores)
        midas_db = MIDAS_DB(os.path.abspath(args.midasdb_dir), args.midasdb_name, num_cores)
        midas_db.fetch_files("pangenome", species_to_analyze)

        if args.prune_centroids:
            centroids_files = {sid: fetch_pruned_centroids(midas_db, sid, args.prune_method, args.prune_cutoff, args.remove_singleton) for sid in species_to_analyze}
        else:
            centroids_files = midas_db.fetch_files("pangenome_centroids", species_to_analyze)

        build_bowtie2_db(bt2_db_dir, bt2_db_name, centroids_files, args.num_cores)
        tsprint("MIDAS2::build_bowtie2db::finish")

        # Align reads to pangenome database
        tsprint("MIDAS2::bowtie2_align::start")
        pangenome_bamfile = sample.get_target_layout("pangenome_bam")
        bowtie2_align(bt2_db_dir, bt2_db_name, pangenome_bamfile, args)
        samtools_index(pangenome_bamfile, args.debug, args.num_cores)
        samtools_idxstats(pangenome_bamfile, args.debug, args.num_cores)
        tsprint("MIDAS2::bowtie2_align::finish")

        if args.alignment_only:
            return

        tsprint("MIDAS2::prepare_species::start")
        prepare_species(species_to_analyze, midas_db)
        tsprint("MIDAS2::prepare_species::finish")

        tsprint("MIDAS2::multiprocessing_map::start")
        # It's important to maintain the order of the bam genes
        fetch_genes_from_bam(f'{pangenome_bamfile}.idxstats', midas_db)

        total_c99_counts = len(readonly_bamgenes)
        number_of_chunks = args.num_cores
        chunk_size = ceil(total_c99_counts / number_of_chunks)

        args_list = []
        chunk_id = 0
        for chunk_start in range(0, total_c99_counts, chunk_size):
            chunk_end = chunk_start + chunk_size
            headerless_sliced_path = sample.get_target_layout("chunk_depth", "", chunk_id)
            args_list.append((chunk_id, chunk_start, chunk_end, pangenome_bamfile, headerless_sliced_path, args.cluster_level))
            chunk_id += 1

        list_of_chunks_depth = multiprocessing_map(compute_pileup_per_chunk, args_list, number_of_chunks)
        dict_of_gene_depth = merge_depth_across_chunks(list_of_chunks_depth, args.cluster_level)
        tsprint("MIDAS2::multiprocessing_map::finish")

        tsprint("MIDAS2::write_species_summary::start")
        species_depth_summary = compute_species_summary(dict_of_gene_depth)
        write_species_summary(species_depth_summary, sample.get_target_layout("genes_summary"))
        tsprint("MIDAS2::write_species_summary::finish")

        if args.remove_bam:
            command(f"rm -f {pangenome_bamfile}", check=False)
            command(f"rm -f {pangenome_bamfile}.bai", check=False)
            command(f"rm -f {pangenome_bamfile}.idxstats", check=False)

        if args.remove_bt2_index:
            sample.remove_dirs(["bt2_indexes_dir"])

        if not args.debug:
            sample.remove_dirs(["tempdir"])

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
