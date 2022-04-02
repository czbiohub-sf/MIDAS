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
from midas2.common.utils import tsprint, InputStream, OutputStream, select_from_tsv, command, multiprocessing_map, multithreading_map, cat_files
from midas2.common.utilities import fetch_genes_are_markers
from midas2.models.midasdb import MIDAS_DB
from midas2.models.sample import Sample
from midas2.models.species import Species, parse_species, load_chunks_cache
from midas2.params.schemas import genes_summary_schema, genes_coverage_schema, format_data, DECIMALS6, genes_chunk_summary_schema, genes_are_markers_schema
from midas2.common.bowtie2 import build_bowtie2_db, bowtie2_align, samtools_index, bowtie2_index_exists, _keep_read
from midas2.params.inputs import MIDASDB_NAMES


DEFAULT_MARKER_DEPTH = 5.0
DEFAULT_MARKER_MEDIAN_DEPTH = 0

DEFAULT_ALN_MAPID = 94.0
DEFAULT_ALN_MAPQ = 2
DEFAULT_ALN_COV = 0.75
DEFAULT_ALN_READQ = 20

DEFAULT_SITE_DEPTH = 2

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
                           help=f"Prebuilt bowtie2 database indexes")
    subparser.add_argument('--prebuilt_bowtie2_species',
                           dest='prebuilt_bowtie2_species',
                           type=str,
                           metavar="CHAR",
                           help=f"List of species used for building the prebuild bowtie2 indexes.")
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
                           help=f"Local MIDAS Database path mirroing S3.")

    subparser.add_argument('--species_list',
                           dest='species_list',
                           type=str,
                           metavar="CHAR",
                           help=f"Comma separated list of species ids")
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
                           help=f"Comman separated correponsding cutoff for filtering species (> XX) ({DEFAULT_MARKER_MEDIAN_DEPTH}, )")

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

    # Site filters
    subparser.add_argument('--site_depth',
                           dest='site_depth',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_SITE_DEPTH,
                           help=f"Minimum number of reads mapped to one gene ({DEFAULT_SITE_DEPTH})")

    # File related
    subparser.add_argument('--chunk_size',
                           dest='chunk_size',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_CHUNK_SIZE,
                           help=f"Number of genomic sites for the temporary chunk file  ({DEFAULT_CHUNK_SIZE})")
    subparser.add_argument('--max_reads',
                           dest='max_reads',
                           type=int,
                           metavar="INT",
                           help=f"Number of reads to use from input file(s).  (All)")
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

    tsprint(f"  MIDAS2::process_chunk_of_genes::{species_id}-{chunk_id}::start compute_coverage_per_chunk")
    ret = compute_coverage_per_chunk(species_id, chunk_id)
    tsprint(f"  MIDAS2::process_chunk_of_genes::{species_id}-{chunk_id}::finish compute_coverage_per_chunk")
    return ret


def compute_coverage_per_chunk(species_id, chunk_id):
    """ Count number of bp mapped to each pan-gene. """

    global global_args
    global semaphore_for_species
    global dict_of_species
    global sample

    try:
        sp = dict_of_species[species_id]

        chunks_of_centroids = load_chunks_cache(sp.chunks_of_centroids_fp)
        dict_of_genes_are_markers, _ = fetch_genes_are_markers(sp.cluster_info_fp)

        chunk_of_genes_length = chunks_of_centroids[chunk_id]
        list_of_centroids_id = sorted(list(chunk_of_genes_length.keys()))

        pangenome_bamfile = sample.get_target_layout("genes_pangenomes_bam")
        headerless_genes_coverage_path = sample.get_target_layout("chunk_coverage", species_id, chunk_id)
        headerless_genes_are_marker = sample.get_target_layout("chunk_genes_are_markers", species_id, chunk_id)

        # Statistics needed to be accmulated within each chunk
        chunk_cov_stats = {
            "species_id": species_id,
            "chunk_id": chunk_id,
            "chunk_genome_size": 0,
            "chunk_num_covered_genes": 0,
            "chunk_coverage": 0,
            "chunk_aligned_reads": 0,
            "chunk_mapped_reads": 0,
        }

        chunk_genes_are_marker = defaultdict(dict)
        chunk_genes_values = defaultdict(list)

        with AlignmentFile(pangenome_bamfile) as bamfile:
            for gene_id in list_of_centroids_id: # compute Unit: Genes
                gene_length = chunk_of_genes_length[gene_id]
                aligned_reads = bamfile.count(gene_id)

                chunk_cov_stats["chunk_genome_size"] += 1
                chunk_cov_stats["chunk_aligned_reads"] += aligned_reads

                mapped_reads = bamfile.count(gene_id, read_callback=keep_read)

                if mapped_reads < global_args.site_depth:
                    continue

                # Compute vertical coverage only for aligned gene region
                gene_covered_bases = 0
                gene_total_depth = 0
                counts = bamfile.count_coverage(gene_id, read_callback=keep_read)
                for within_chunk_index in range(0, gene_length):
                    site_depth = sum([counts[nt][within_chunk_index] for nt in range(4)])
                    gene_total_depth += site_depth
                    if site_depth > 0:
                        gene_covered_bases += 1

                gene_fraction_covered = gene_covered_bases / gene_length
                gene_coverage = gene_total_depth / gene_covered_bases if gene_covered_bases > 0 else 0

                if gene_coverage == 0: # Sparse by default.
                    continue

                chunk_cov_stats["chunk_num_covered_genes"] += 1
                chunk_cov_stats["chunk_mapped_reads"] += mapped_reads
                chunk_cov_stats["chunk_coverage"] += gene_coverage

                chunk_genes_values[gene_id] = [gene_id, gene_length, aligned_reads, mapped_reads, gene_coverage, gene_fraction_covered, 0.0]

                # Collect gene coverage for all the centroids (genes) that are markers
                if gene_id in dict_of_genes_are_markers:
                    chunk_genes_are_marker[gene_id] = dict_of_genes_are_markers[gene_id]
                    chunk_genes_are_marker[gene_id]["gene_coverage"] = gene_coverage

        # Write chunk gene coverage results to file
        with OutputStream(headerless_genes_coverage_path) as stream:
            for vals in chunk_genes_values.values():
                stream.write("\t".join(map(format_data, vals, repeat(DECIMALS6, len(vals)))) + "\n")
        # Write current chunk's genes_that_are_markers to file
        with OutputStream(headerless_genes_are_marker) as stream2:
            for rec in chunk_genes_are_marker.values():
                vals = rec.values()
                stream2.write("\t".join(map(format_data, vals, repeat(DECIMALS6, len(vals)))) + "\n")

        return chunk_cov_stats
    finally:
        semaphore_for_species[species_id].release()


def compute_scg_coverage_across_chunks(species_id):
    """ Extract gene depth for mapped centroids_99 that are markers and Compute the median marker coverage """

    global dict_of_species
    sp = dict_of_species[species_id]

    number_of_chunks = sp.num_of_genes_chunks
    _, list_of_markers = fetch_genes_are_markers(sp.cluster_info_fp)

    # Include ALL the marker centroids, not only the covered one
    markers_coverage = dict(zip(list_of_markers, [0.0]*len(list_of_markers)))

    list_of_chunk_files = [sample.get_target_layout("chunk_genes_are_markers", species_id, chunk_id) for chunk_id in range(0, number_of_chunks)]
    for chunk_file in list_of_chunk_files:
        with InputStream(chunk_file) as stream:
            for r in select_from_tsv(stream, schema=genes_are_markers_schema, result_structure=dict):
                if r["marker_id"] not in markers_coverage:
                    markers_coverage[r["marker_id"]] = 0.0
                markers_coverage[r["marker_id"]] += r["gene_coverage"]

    median_marker_coverage = np.median(list(map(lambda x: float(format(x, DECIMALS6)), markers_coverage.values())))

    return median_marker_coverage


def update_chunk_coverage(my_args):

    chunk_coverage_path, median_marker_coverage = my_args

    c_copynum = list(genes_coverage_schema.keys()).index("copy_number")
    c_coverage = list(genes_coverage_schema.keys()).index("mean_coverage")

    add_cn_to_write = []
    with InputStream(chunk_coverage_path) as stream:
        for line in stream:
            vals = line.rstrip("\n").split("\t")
            # Infer gene copy counts
            vals[c_copynum] = float(vals[c_coverage]) / median_marker_coverage
            add_cn_to_write.append(vals)

    with OutputStream(chunk_coverage_path) as stream:
        for vals in add_cn_to_write:
            stream.write("\t".join(map(format_data, vals, repeat(DECIMALS6, len(vals)))) + "\n")


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

    # Overwrite the chunk_gene_coverage file with updated copy_number
    list_of_chunks_coverage = [sample.get_target_layout("chunk_coverage", species_id, chunk_id) for chunk_id in range(0, number_of_chunks)]
    species_gene_coverage_path = sample.get_target_layout("genes_coverage", species_id)

    tsprint(f"      MIDAS2::merge_chunks_per_species::{species_id}::start update_chunk_coverage")
    if median_marker_coverage > 0:
        args = [(chunk_file, median_marker_coverage) for chunk_file in list_of_chunks_coverage]
        multithreading_map(update_chunk_coverage, args, 4)
    tsprint(f"      MIDAS2::merge_chunks_per_species::{species_id}::finish update_chunk_coverage")

    # Merge chunks results into per species genes coverage file
    with OutputStream(species_gene_coverage_path) as stream:
        stream.write('\t'.join(genes_coverage_schema.keys()) + '\n')
    cat_files(list_of_chunks_coverage, species_gene_coverage_path, 10)

    if not global_args.debug:
        tsprint(f"Deleting temporary sliced coverage files for {species_id}.")
        temp_species_dir = os.path.dirname(list_of_chunks_coverage[0])
        command(f"rm -rf {temp_species_dir}", quiet=True)

    return {"species_id": species_id, "chunk_id": -1, "median_marker_coverage": median_marker_coverage}


def write_species_coverage_summary(chunks_gene_coverage, genes_stats_path):

    # First pass
    species_coverage_summary = defaultdict(dict)
    for rec in chunks_gene_coverage:
        species_id = rec["species_id"]

        # for the merge task, we return the marker genes coverage
        if rec["chunk_id"] == -1:
            species_coverage_summary[species_id]["median_marker_coverage"] = rec["median_marker_coverage"]
            continue

        if species_id not in species_coverage_summary:
            species_coverage_summary[species_id] = {
                "species_id": species_id,
                "pangenome_size": 0,
                "num_covered_genes": 0,
                "total_coverage": 0,
                "aligned_reads": 0,
                "mapped_reads": 0,
                "median_marker_coverage": 0,
            }

        species_coverage_summary[species_id]["pangenome_size"] += rec["chunk_genome_size"]
        species_coverage_summary[species_id]["num_covered_genes"] += rec["chunk_num_covered_genes"]
        species_coverage_summary[species_id]["aligned_reads"] += rec["chunk_aligned_reads"]
        species_coverage_summary[species_id]["mapped_reads"] += rec["chunk_mapped_reads"]
        species_coverage_summary[species_id]["total_coverage"] += rec["chunk_coverage"]

    # Second pass: need to loop over all the chunks to calculate the average read coverage
    with OutputStream(genes_stats_path) as stream:
        stream.write("\t".join(genes_summary_schema.keys()) + "\n")
        for rec in species_coverage_summary.values():
            fraction_covered = rec["num_covered_genes"] / rec["pangenome_size"]
            mean_coverage = rec["total_coverage"] / rec["num_covered_genes"] if rec["num_covered_genes"] > 0 else 0

            vals = [rec["species_id"], rec["pangenome_size"], rec["num_covered_genes"], \
                    fraction_covered, mean_coverage, rec["aligned_reads"], rec["mapped_reads"], rec["median_marker_coverage"]]
            stream.write("\t".join(map(format_data, vals)) + "\n")


def write_chunk_coverage_summary(chunks_gene_coverage, outfile):
    with OutputStream(outfile) as stream:
        stream.write("\t".join(genes_chunk_summary_schema.keys()) + "\n")
        for rec in chunks_gene_coverage:
            if rec["chunk_id"] == -1:
                continue
            stream.write("\t".join(map(format_data, rec.values())) + "\n")


def run_genes(args):

    try:
        global global_args
        global_args = args

        global sample
        sample = Sample(args.sample_name, args.midas_outdir, "genes")
        sample.create_dirs(["outdir", "tempdir"], args.debug)

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

        arguments_list = design_chunks(species_ids_of_interest, midas_db, args.chunk_size)
        tsprint(f"MIDAS2::design_chunks::finish")

        # Build Bowtie indexes for species in the restricted species profile
        centroids_files = midas_db.fetch_files("pangenome_centroids", species_ids_of_interest)
        midas_db.fetch_files("pangenome_cluster_info", species_ids_of_interest)
        tsprint(f"MIDAS2::build_bowtie2db::start")
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
        write_chunk_coverage_summary(chunks_gene_coverage, sample.get_target_layout("genes_chunk_summary"))
        tsprint(f"MIDAS2::write_species_coverage_summary::finish")

    except Exception as error:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error.  Specify --debug flag to keep.")
            sample.remove_dirs(["outdir", "tempdir"])
            if not args.prebuilt_bowtie2_indexes:
                sample.remove_dirs(["bt2_indexes_dir"])
        raise error


@register_args
def main(args):
    tsprint(f"Pangenome copy number profiling in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    run_genes(args)
