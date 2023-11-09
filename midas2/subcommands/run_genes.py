#!/usr/bin/env python3
import json
import os
import multiprocessing
from operator import itemgetter

from collections import defaultdict
from itertools import repeat
import numpy as np
from pysam import AlignmentFile  # pylint: disable=no-name-in-module

from midas2.common.argparser import add_subcommand
from midas2.common.utils import tsprint, InputStream, OutputStream, select_from_tsv, command, multiprocessing_map, multithreading_map, cat_files, args_string
from midas2.common.utilities import fetch_c99s_are_markers
from midas2.models.midasdb import MIDAS_DB
from midas2.models.sample import Sample
from midas2.models.species import Species, parse_species, load_chunks_cache, scan_cluster_info
from midas2.params.schemas import genes_summary_schema, genes_coverage_schema, format_data, DECIMALS6, genes_are_markers_schema, species_marker_profile_schema
from midas2.common.bowtie2 import build_bowtie2_db, bowtie2_align, samtools_index, bowtie2_index_exists, _keep_read
from midas2.params.inputs import MIDASDB_NAMES


DEFAULT_MARKER_DEPTH = 5.0
DEFAULT_MARKER_MEDIAN_DEPTH = 2.0

DEFAULT_ALN_MAPID = 94.0
DEFAULT_ALN_MAPQ = 0
DEFAULT_ALN_COV = 0.75
DEFAULT_ALN_READQ = 20

DEFAULT_READ_DEPTH = 0

DEFAULT_CHUNK_SIZE = 50000
DEFAULT_MAX_FRAGLEN = 5000
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
                           help=f"Path to prebuilt pan-genome bowtie2 database indexes")
    subparser.add_argument('--prebuilt_bowtie2_species',
                           dest='prebuilt_bowtie2_species',
                           type=str,
                           metavar="CHAR",
                           help=f"Path to list of species in the prebuild bowtie2 indexes.")
    subparser.add_argument('--midasdb_name',
                           dest='midasdb_name',
                           type=str,
                           default="uhgg",
                           choices=MIDASDB_NAMES,
                           help=f"MIDAS Database name.")
    subparser.add_argument('--midasdb_dir',
                           dest='midasdb_dir',
                           type=str,
                           default="midasdb",
                           help=f"Path to local MIDAS Database.")

    subparser.add_argument('--species_list',
                           dest='species_list',
                           type=str,
                           metavar="CHAR",
                           help=f"Comma separated list of species ids OR path to list of species TXT.")
    subparser.add_argument('--select_by',
                           dest='select_by',
                           type=str,
                           default="median_marker_coverage",
                           help=f"Comma separated columns from species_profile to filter species.")
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
                           help=f"Maximum fragment length for paired reads.")
    subparser.add_argument('--max_reads',
                           dest='max_reads',
                           type=int,
                           metavar="INT",
                           help=f"Number of reads to use from input file(s) for read alignment.  (All)")

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
    subparser.add_argument('--chunk_size',
                           dest='chunk_size',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_CHUNK_SIZE,
                           help=f"Number of genomic sites for the temporary chunk file  ({DEFAULT_CHUNK_SIZE})")
    subparser.add_argument('--num_cores',
                           dest='num_cores',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_NUM_CORES,
                           help=f"Number of physical cores to use ({DEFAULT_NUM_CORES})")
    return main_func


def keep_read(aln):
    global global_args
    args = global_args
    return _keep_read(aln, args.aln_mapid, args.aln_readq, args.aln_mapq, args.aln_cov)


def design_chunks_per_species(pargs):
    """ Fetch given species's pre-computed gene_id - centroid_99 - marker_id mapping information """
    sp, midas_db, chunk_size = pargs
    return sp.compute_gene_chunks(midas_db, chunk_size)


def design_chunks(species_ids_of_interest, midas_db, chunk_size):
    """ Each chunk_of_genes is indexed by (species_id, chunk_id) """

    global dict_of_species
    global semaphore_for_species

    semaphore_for_species = dict()
    dict_of_species = {species_id: Species(species_id) for species_id in species_ids_of_interest}

    num_cores = min(midas_db.num_cores, 10)
    flags = multithreading_map(design_chunks_per_species, [(sp, midas_db, chunk_size) for sp in dict_of_species.values()], num_cores) #<---
    assert all(flags)

    # Sort species by ascending num_of_genes_chunks
    sorted_tuples_of_species = sorted(((sp.id, sp.num_of_genes_chunks) for sp in dict_of_species.values()), key=itemgetter(1))

    arguments_list = []
    for species_id, num_of_genes_chunks in sorted_tuples_of_species:
        sp = dict_of_species[species_id]

        for chunk_id in range(0, num_of_genes_chunks):
            arguments_list.append((species_id, chunk_id))

        semaphore_for_species[species_id] = multiprocessing.Semaphore(num_of_genes_chunks)
        for _ in range(num_of_genes_chunks):
            semaphore_for_species[species_id].acquire()

    tsprint("================= Total number of compute chunks: " + str(len(arguments_list)))

    # Submit the merge tasks to the end of the queue.
    for species_id, sp in dict_of_species.items():
        arguments_list.append((species_id, -1, sp.num_of_genes_chunks))
    return arguments_list


def process_chunk_of_genes(packed_args):
    """ Control flow of compute coverage of pangenome per species and write results """

    species_id, chunk_id = packed_args[:2]

    if chunk_id == -1:
        global semaphore_for_species
        num_of_genes_chunks = packed_args[2]

        tsprint(f"  MIDAS2::process_chunk_of_genes::{species_id}-{chunk_id}::wait merge_chunks_per_species")
        for _ in range(num_of_genes_chunks):
            semaphore_for_species[species_id].acquire()
        tsprint(f"  MIDAS2::process_chunk_of_genes::{species_id}-{chunk_id}::start merge_chunks_per_species")
        ret = merge_chunks_per_species(species_id)
        tsprint(f"  MIDAS2::process_chunk_of_genes::{species_id}-{chunk_id}::finish merge_chunks_per_species")

        return ret

    tsprint(f"  MIDAS2::process_chunk_of_genes::{species_id}-{chunk_id}::start collect_reads_per_chunk")
    ret = collect_reads_per_chunk(species_id, chunk_id)
    tsprint(f"  MIDAS2::process_chunk_of_genes::{species_id}-{chunk_id}::finish collect_reads_per_chunk")
    return ret


def collect_reads_per_chunk(species_id, chunk_id):
    """ Collect total number of read depths over each covered position for centroids_99 """

    global global_args
    global semaphore_for_species
    global dict_of_species
    global sample

    try:
        sp = dict_of_species[species_id]

        chunks_of_centroids = load_chunks_cache(sp.chunks_of_centroids_fp)
        dict_of_c99s_are_markers, _ = fetch_c99s_are_markers(sp.cluster_info_fp)
        cluster_info = scan_cluster_info(sp.cluster_info_fp)

        chunk_of_c99_length = chunks_of_centroids[chunk_id]
        list_of_c99_ids = sorted(list(chunk_of_c99_length.keys()))

        pangenome_bamfile = sample.get_target_layout("genes_pangenomes_bam")
        headerless_c95_are_marker = sample.get_target_layout("chunk_c95_are_markers", species_id, chunk_id)
        headerless_c95_coverage_path = sample.get_target_layout("chunk_c95_coverage", species_id, chunk_id)

        chunk_c95_values = defaultdict(dict)
        chunk_c95_are_marker = defaultdict(dict)
        with AlignmentFile(pangenome_bamfile) as bamfile:
            # Competitive alignment is done on centroid_99 level, while coverage are computed on centroid_95 level.
            for c99_id in list_of_c99_ids:
                c95_id = cluster_info[c99_id]["centroid_95"]

                # centroids_95 is subset of centroids_99
                c99_length = cluster_info[c99_id]["centroid_99_length"]
                c95_length = cluster_info[c95_id]["centroid_99_length"]

                aligned_reads = bamfile.count(c99_id)
                mapped_reads = bamfile.count(c99_id, read_callback=keep_read)

                # Compute total per-position depth for aligned gene region
                c99_covered_bases = 0
                c99_total_depth = 0
                counts = bamfile.count_coverage(c99_id, read_callback=keep_read)
                for within_chunk_index in range(0, c99_length):
                    # Per-position depth
                    gene_depth = sum([counts[nt][within_chunk_index] for nt in range(4)])
                    c99_total_depth += gene_depth
                    if gene_depth > 0:
                        c99_covered_bases += 1

                if c99_total_depth == 0: # Sparse by default.
                    continue

                if c95_id not in chunk_c95_values:
                    chunk_c95_values[c95_id]["c95_id"] = c95_id
                    chunk_c95_values[c95_id]["c95_length"] = c95_length
                    chunk_c95_values[c95_id]["aligned_reads"] = aligned_reads
                    chunk_c95_values[c95_id]["mapped_reads"] = mapped_reads
                    chunk_c95_values[c95_id]["read_depth"] = c99_total_depth
                    chunk_c95_values[c95_id]["mean_coverage"] = 0.0
                    chunk_c95_values[c95_id]["copy_number"] = 0.0
                else:
                    chunk_c95_values[c95_id]["aligned_reads"] += aligned_reads
                    chunk_c95_values[c95_id]["mapped_reads"] += mapped_reads
                    chunk_c95_values[c95_id]["read_depth"] += c99_total_depth

                # Collect gene coverage for all centroids_95 that are markers, or member of the centroids_99 is marker
                if c95_id in dict_of_c99s_are_markers and c95_id not in chunk_c95_are_marker:
                    chunk_c95_are_marker[c95_id]["centroid_95"] = c95_id
                    chunk_c95_are_marker[c95_id]["marker_id"] = dict_of_c99s_are_markers[c95_id]["marker_id"]
                    chunk_c95_are_marker[c95_id]["marker_length"] = c95_length
                    chunk_c95_are_marker[c95_id]["read_depth"] = c99_total_depth
                elif c95_id in chunk_c95_are_marker:
                    chunk_c95_are_marker[c95_id]["read_depth"] += c99_total_depth

        # Write chunk centroids_95 that has at least one marker member to file
        with OutputStream(headerless_c95_are_marker) as stream:
            for rec in chunk_c95_are_marker.values():
                vals = rec.values()
                stream.write("\t".join(map(format_data, vals, repeat(DECIMALS6, len(vals)))) + "\n")

        # Write chunk centroid_95 coverage to file
        with OutputStream(headerless_c95_coverage_path) as stream:
            for rec in chunk_c95_values.values():
                vals = rec.values()
                stream.write("\t".join(map(format_data, vals, repeat(DECIMALS6, len(vals)))) + "\n")

        return {"species_id": species_id, "chunk_id": chunk_id}
    finally:
        semaphore_for_species[species_id].release()


def compute_scg_coverage_across_chunks(species_id):
    """ Extract gene depth for mapped centroids_99 that are markers and Compute the median marker coverage """

    global dict_of_species
    sp = dict_of_species[species_id]

    number_of_chunks = sp.num_of_genes_chunks
    _, list_of_markers = fetch_c99s_are_markers(sp.cluster_info_fp)

    # Gene flow: ALL; Species flow: only covered markers.
    list_of_chunks_marker = [sample.get_target_layout("chunk_c95_are_markers", species_id, chunk_id) for chunk_id in range(0, number_of_chunks)]
    markers_coverage = dict(zip(list_of_markers, [0.0]*len(list_of_markers)))
    for chunk_file in list_of_chunks_marker:
        with InputStream(chunk_file) as stream:
            for r in select_from_tsv(stream, schema=genes_are_markers_schema, result_structure=dict):
                if r["marker_id"] not in markers_coverage:
                    markers_coverage[r["marker_id"]] = 0.0
                markers_coverage[r["marker_id"]] += r["read_depth"] / r["marker_length"] # alias for c95_length

    median_marker_coverage = np.median(list(map(lambda x: float(format(x, DECIMALS6)), markers_coverage.values())))
    return median_marker_coverage


def compute_coverage_across_chunks(list_of_chunks_coverage, median_marker_coverage):
    # Second round of across chunks compute: same centroid_95 would be present in differnet chunks.
    # Therefore, have to be scanned sequentially.
    c95_values = defaultdict(dict)
    for chunk_file in list_of_chunks_coverage:
        with InputStream(chunk_file) as stream:
            for r in select_from_tsv(stream, schema=genes_coverage_schema, result_structure=dict):
                if r["gene_id"] not in c95_values:
                    c95_values[r["gene_id"]] = r
                else:
                    c95_values[r["gene_id"]]["aligned_reads"] += r["aligned_reads"]
                    c95_values[r["gene_id"]]["mapped_reads"] += r["mapped_reads"]
                    c95_values[r["gene_id"]]["read_depth"] += r["read_depth"]

    # Compute coverage and copy number
    for c95_id, c95_cov in c95_values.items():
        c95_cov["mean_coverage"] = float(c95_cov["read_depth"] / c95_cov["gene_length"])
        if median_marker_coverage > 0:
            c95_cov["copy_number"] = float(c95_cov["mean_coverage"] / median_marker_coverage)
    return c95_values


def merge_chunks_per_species(species_id):
    """ Compute coverage of pangenome for given species_id and write results to disk """

    global sample
    global dict_of_species
    global global_args

    sp = dict_of_species[species_id]
    number_of_chunks = sp.num_of_genes_chunks

    # Compute the median read coverage for all mapped marker genes for given species
    tsprint(f"      MIDAS2::merge_chunks_per_species::{species_id}::start compute_scg_coverage_across_chunks")
    median_marker_coverage = compute_scg_coverage_across_chunks(species_id)
    tsprint(f"      MIDAS2::merge_chunks_per_species::{species_id}::finish compute_scg_coverage_across_chunks")

    # Compute centroid_95 coverage and copy number
    tsprint(f"      MIDAS2::merge_chunks_per_species::{species_id}::start compute_coverage_across_chunks")
    list_of_chunks_coverage = [sample.get_target_layout("chunk_c95_coverage", species_id, chunk_id) for chunk_id in range(0, number_of_chunks)]
    c95_vals = compute_coverage_across_chunks(list_of_chunks_coverage, median_marker_coverage)
    tsprint(f"      MIDAS2::merge_chunks_per_species::{species_id}::finish compute_coverage_across_chunks")

    species_cov_stats = {
        "species_id": species_id,
        "chunk_id": -1,
        "covered_genes": 0,
        "mean_coverage": 0,
        "aligned_reads": 0,
        "mapped_reads": 0,
        "marker_coverage": median_marker_coverage
    }

    # Write centroid_95 coverage to file.
    species_gene_cov_path = sample.get_target_layout("genes_coverage", species_id)
    with OutputStream(species_gene_cov_path) as stream:
        stream.write('\t'.join(genes_coverage_schema.keys()) + '\n')
        for rec in c95_vals.values():
            if rec["mapped_reads"] > global_args.read_depth:
                species_cov_stats["covered_genes"] += 1
                species_cov_stats["mean_coverage"] += rec["mean_coverage"]
                species_cov_stats["aligned_reads"] += rec["aligned_reads"]
                species_cov_stats["mapped_reads"] += rec["mapped_reads"]
                stream.write("\t".join(map(format_data, rec.values(), repeat(DECIMALS6, len(rec)))) + "\n")

    # Remove chunk temporary file TODO remove this
    tsprint(f"Deleting temporary chunk coverage files for {species_id}.")
    temp_species_dir = os.path.dirname(list_of_chunks_coverage[0])
    command(f"rm -rf {temp_species_dir}", quiet=True)

    return species_cov_stats


def write_species_coverage_summary(chunks_gene_coverage, genes_coverage_path):
    genes_summary = defaultdict(dict)
    for rec in chunks_gene_coverage:
        if rec["chunk_id"] == -1:
            spid = rec["species_id"]
            sp = dict_of_species[spid]
            cluster_info = scan_cluster_info(sp.cluster_info_fp)

            # Compute the number of unique centroid_95 across all chunks
            c95_dict = dict()
            for ci in cluster_info.values():
                c95 = ci["centroid_95"]
                if c95 not in c95_dict:
                    c95_dict[c95] = 1
                else:
                    c95_dict[c95] += 1
            pangenome_size = len(c95_dict)

            genes_summary[spid]["species_id"] = spid
            genes_summary[spid]["pangenome_size"] = pangenome_size
            genes_summary[spid]["covered_genes"] = rec["covered_genes"]
            genes_summary[spid]["fraction_covered"] = rec["covered_genes"] / pangenome_size
            genes_summary[spid]["aligned_reads"] = rec["aligned_reads"]
            genes_summary[spid]["mapped_reads"] = rec["mapped_reads"]
            genes_summary[spid]["mean_coverage"] = rec["mean_coverage"] / rec["covered_genes"] if rec["covered_genes"] > 0 else 0
            genes_summary[spid]["marker_coverage"] = rec["marker_coverage"]
            continue
    
    with OutputStream(genes_coverage_path) as stream:
        stream.write("\t".join(genes_summary_schema.keys()) + "\n")
        for _ in genes_summary.values():
            #vals = [rec["species_id"], rec["pangenome_size"], rec["num_covered_genes"], \
            #        fraction_covered, mean_coverage, rec["aligned_reads"], rec["mapped_reads"], rec["median_marker_coverage"]]
            stream.write("\t".join(map(format_data, _.values())) + "\n")


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
            assert (args.prebuilt_bowtie2_species and os.path.exists(args.prebuilt_bowtie2_species)), f"Require list of speices used to build the provided Bowtie2 indexes."

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

        species_ids_of_interest = species_list if no_filter else sample.select_species(args, species_list)
        species_counts = len(species_ids_of_interest)

        assert species_counts > 0, f"No (specified) species pass the marker_depth filter, please adjust the marker_depth or species_list"
        tsprint(len(species_ids_of_interest))

        tsprint(f"MIDAS2::design_chunks::start")
        num_cores_download = min(species_counts, args.num_cores)
        midas_db = MIDAS_DB(os.path.abspath(args.midasdb_dir), args.midasdb_name, num_cores_download)
        midas_db.fetch_files("pangenome", species_ids_of_interest)
        arguments_list = design_chunks(species_ids_of_interest, midas_db, args.chunk_size)
        tsprint(f"MIDAS2::design_chunks::finish")

        # Build Bowtie indexes for species in the restricted species profile
        tsprint(f"MIDAS2::build_bowtie2db::start")
        centroids_files = midas_db.fetch_files("pangenome_centroids", species_ids_of_interest)
        build_bowtie2_db(bt2_db_dir, bt2_db_name, centroids_files, args.num_cores)
        tsprint(f"MIDAS2::build_bowtie2db::finish")

        # Align reads to pangenome database
        tsprint(f"MIDAS2::bowtie2_align::start")
        sample.create_species_subdirs(species_ids_of_interest, "temp", args.debug)
        pangenome_bamfile = sample.get_target_layout("genes_pangenomes_bam")
        bowtie2_align(bt2_db_dir, bt2_db_name, pangenome_bamfile, args)
        samtools_index(pangenome_bamfile, args.debug, args.num_cores)
        tsprint(f"MIDAS2::bowtie2_align::finish")

        tsprint(f"MIDAS2::multiprocessing_map::start")
        total_chunk_counts = sum((sp.num_of_genes_chunks for sp in dict_of_species.values()))
        num_cores = min(args.num_cores, total_chunk_counts)
        tsprint(f"The number of cores will be used to compute coverage is {num_cores}")
        chunks_gene_coverage = multiprocessing_map(process_chunk_of_genes, arguments_list, num_cores)
        tsprint(f"MIDAS2::multiprocessing_map::finish")

        tsprint(f"MIDAS2::write_species_coverage_summary::start")
        write_species_coverage_summary(chunks_gene_coverage, sample.get_target_layout("genes_summary"))
        tsprint(f"MIDAS2::write_species_coverage_summary::finish")

        if not args.debug:
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
