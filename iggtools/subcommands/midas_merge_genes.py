#!/usr/bin/env python3
import os
import json
from collections import defaultdict
import multiprocessing
from itertools import repeat

from iggtools.models.samplepool import SamplePool
from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, InputStream, OutputStream, select_from_tsv, multiprocessing_map, num_physical_cores, cat_files
from iggtools.models.midasdb import MIDAS_DB
from iggtools.params.schemas import genes_info_schema, genes_coverage_schema, format_data, fetch_default_genome_depth, DECIMALS6


DEFAULT_GENOME_DEPTH = fetch_default_genome_depth("genes")
DEFAULT_SAMPLE_COUNTS = 1
DEFAULT_CLUSTER_ID = '95'
DEFAULT_MIN_COPY = 0.35
DEFAULT_CHUNK_SIZE = 50000


def register_args(main_func):
    subparser = add_subcommand('midas_merge_genes', main_func, help='metagenomic pan-genome profiling')

    subparser.add_argument('midas_outdir',
                           type=str,
                           help="""Path to directory to store results.  Subdirectory will be created for each species.""")
    subparser.add_argument('--samples_list',
                           dest='samples_list',
                           type=str,
                           required=True,
                           help=f"TSV file mapping sample name to midas_run_species.py output directories")
    subparser.add_argument('--chunk_size',
                           dest='chunk_size',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_CHUNK_SIZE,
                           help=f"Number of genomic sites for the temporary chunk file  ({DEFAULT_CHUNK_SIZE})")

    # Species and sample filters
    subparser.add_argument('--species_list',
                           dest='species_list',
                           type=str,
                           metavar="CHAR",
                           help=f"Comma separated list of species ids")
    subparser.add_argument('--genome_depth',
                           dest='genome_depth',
                           type=float,
                           metavar="FLOAT",
                           default=DEFAULT_GENOME_DEPTH,
                           help=f"Minimum average read depth per sample ({DEFAULT_GENOME_DEPTH})")
    subparser.add_argument('--sample_counts',
                           dest='sample_counts', #min_samples
                           type=int,
                           metavar="INT",
                           default=DEFAULT_SAMPLE_COUNTS,
                           help=f"select species with >= MIN_SAMPLES ({DEFAULT_SAMPLE_COUNTS})")

    subparser.add_argument('--midas_db',
                           dest='midas_db',
                           type=str,
                           metavar="CHAR",
                           help=f"local MIDAS DB which mirrors the s3 IGG db")
    subparser.add_argument('--num_cores',
                           dest='num_cores',
                           type=int,
                           metavar="INT",
                           default=num_physical_cores,
                           help=f"Number of physical cores to use ({num_physical_cores})")

    # Presence/Absence
    subparser.add_argument('--min_copy',
                           dest='min_copy',
                           type=float,
                           metavar="FLOAT",
                           default=DEFAULT_MIN_COPY,
                           help=f"Genes >= MIN_COPY are classified as present ({DEFAULT_MIN_COPY})")
    subparser.add_argument('--cluster_pid',
                           dest='cluster_pid',
                           type=str,
                           default=DEFAULT_CLUSTER_ID,
                           choices=['75', '80', '85', '90', '95', '99'],
                           help=f"CLUSTER_PID allows you to quantify gene content for any of these sets of gene clusters ({DEFAULT_CLUSTER_ID})")

    return main_func


def design_chunks(midas_db, chunk_size):
    """ Each chunk_of_genes is indexed by (species_id, chunk_id) """

    global dict_of_species
    global semaphore_for_species
    semaphore_for_species = dict()

    # Design chunks of genes structure per species
    #assert all([design_chunks_per_species((sp, midas_db, chunk_size)) for sp in dict_of_species.values()])
    # sorted_tuples_of_species = sorted(((sp.id, sp.num_of_genes_chunks) for sp in dict_of_species.values()), key=itemgetter(1))
    #for species_id, num_of_genes_chunks in sorted_tuples_of_species:

    arguments_list = []
    for species_id, sp in dict_of_species.items():
        species_id = sp.id

        assert sp.design_genes_chunks(midas_db, chunk_size)
        num_of_chunks = sp.num_of_genes_chunks

        for chunk_id in range(0, num_of_chunks):
            arguments_list.append((species_id, chunk_id))

        semaphore_for_species[species_id] = multiprocessing.Semaphore(num_of_chunks)
        for _ in range(num_of_chunks):
            semaphore_for_species[species_id].acquire()

    #for species_id in dict_of_species.keys():
    #    arguments_list.append((species_id, -1))

    return arguments_list


def process_one_chunk_of_genes(packed_args):

    species_id, chunk_id = packed_args

    if chunk_id == -1:
        global semaphore_for_species
        global dict_of_species

        sp = dict_of_species[species_id]

        tsprint(f"  CZ::process_one_chunk_of_genes::{species_id}--1::wait merge_all_chunks_per_species")
        for _ in range(sp.num_of_genes_chunks):
            semaphore_for_species[species_id].acquire()
        tsprint(f"  CZ::process_one_chunk_of_genes::{species_id}--1::start merge_all_chunks_per_species")
        merge_all_chunks_per_species(species_id)
        tsprint(f"  CZ::process_one_chunk_of_genes::{species_id}--1::finish merge_all_chunks_per_species")
        return "worked"

    # Write to chunk files
    tsprint(f"  CZ::process_one_chunk_of_genes::{species_id}-{chunk_id}::start pool_across_samples_per_chunk")
    pool_across_samples_per_chunk(species_id, chunk_id)
    tsprint(f"  CZ::process_one_chunk_of_genes::{species_id}-{chunk_id}::finish pool_across_samples_per_chunk")

    return "worked"


def pool_across_samples_per_chunk(species_id, chunk_id):

    global semaphore_for_species
    global pool_of_samples
    global dict_of_species

    global global_args

    try:
        sp = dict_of_species[species_id]
        total_samples_count = sp.samples_count

        # First pass: accumulate the matrix sample by sample
        list_of_genes_id = set([sp.clusters_map[g_id] for g_id in sp.chunks_of_centroids[chunk_id].keys()])
        accumulator = {
                        "copynum": {g_id:[0.0] * total_samples_count for g_id in list_of_genes_id},
                        "depth":   {g_id:[0.0] * total_samples_count for g_id in list_of_genes_id},
                        "reads":   {g_id:[0.0] * total_samples_count for g_id in list_of_genes_id},
                        "presabs": {g_id:[0] * total_samples_count for g_id in list_of_genes_id}
                      }

        tsprint(f"    CZ::process_one_chunk_of_genes::{species_id}-{chunk_id}::start collect_sample_by_sample")
        for sample_index, sample in enumerate(sp.list_of_samples):
            species_coverage_path = sample.get_target_layout("genes_coverage", species_id)
            chunk_coverage_path = sample.get_target_layout("chunk_coverage", species_id, chunk_id)

            genes_coverage_path = chunk_coverage_path if os.path.exists(chunk_coverage_path) else species_coverage_path
            has_header = not os.path.exists(chunk_coverage_path)

            tsprint(f"    CZ::process_one_chunk_of_genes::{species_id}-{chunk_id}::start collect_one_sample")
            my_args = (species_id, chunk_id, sample_index, genes_coverage_path, has_header)
            collect(accumulator, my_args)
            tsprint(f"    CZ::process_one_chunk_of_genes::{species_id}-{chunk_id}::finish collect_one_sample")
        tsprint(f"    CZ::process_one_chunk_of_genes::{species_id}-{chunk_id}::finish collect_sample_by_sample")

        # Second pass: infer presence absence based on copy number
        for gene_id, copynum in accumulator["copynum"].items():
            accumulator["presabs"][gene_id] = [1 if cn >= global_args.min_copy else 0 for cn in copynum]

        # Write to chunk file
        tsprint(f"    CZ::process_one_chunk_of_genes::{species_id}-{chunk_id}::start write_pooled_genes_per_chunk")
        assert write_pooled_genes_per_chunk(accumulator, species_id, chunk_id)
        tsprint(f"    CZ::process_one_chunk_of_genes::{species_id}-{chunk_id}::finish write_pooled_genes_per_chunk")
    finally:
        semaphore_for_species[species_id].release() # no deadlock


def collect(accumulator, my_args):

    global dict_of_species

    species_id, chunk_id, sample_index, genes_coverage_path, has_header = my_args

    sp = dict_of_species[species_id]
    clusters_map = sp.clusters_map
    list_of_centroids_id = sorted(list(sp.chunks_of_centroids[chunk_id].keys()))

    with InputStream(genes_coverage_path) as stream:
        if has_header:
            stream.readline()

        for r in select_from_tsv(stream, schema=genes_coverage_schema, result_structure=dict):
            if r["gene_id"] not in list_of_centroids_id:
                continue

            gene_id = clusters_map[r["gene_id"]]
            accumulator["copynum"][gene_id][sample_index] += r["copy_number"]
            accumulator["depth"][gene_id][sample_index] += r["total_depth"]
            accumulator["reads"][gene_id][sample_index] += r["mapped_reads"]


def write_pooled_genes_per_chunk(accumulator, species_id, chunk_id):
    global pool_of_samples
    for file_type in list(genes_info_schema.keys()):
        outfile = pool_of_samples.get_target_layout(f"genes_{file_type}_by_chunk", species_id, chunk_id)
        with OutputStream(outfile) as stream:
            for gene_id, gene_vals in accumulator[file_type].items():
                stream.write(f"{gene_id}\t" + "\t".join(map(format_data, gene_vals, repeat(DECIMALS6, len(gene_vals)))) + "\n")
    return True


def merge_all_chunks_per_species(species_id):
    global global_args
    global dict_of_species
    global pool_of_samples

    sp = dict_of_species[species_id]
    number_of_chunks = sp.num_of_genes_chunks
    samples_names = dict_of_species[species_id].fetch_samples_names()

    for file_type in list(genes_info_schema.keys()):
        list_of_chunks_path = [pool_of_samples.get_target_layout(f"genes_{file_type}_by_chunk", species_id, chunk_id) for chunk_id in range(0, number_of_chunks)]
        species_path = pool_of_samples.get_target_layout(f"genes_{file_type}", species_id)

        with OutputStream(species_path) as stream:
            stream.write("gene_id\t" + "\t".join(samples_names) + "\n")
        cat_files(list_of_chunks_path, species_path, 10)

        #if not global_args.debug:
        #    for s_file in list_of_chunks_path:
        #        command(f"rm -rf {s_file}", quiet=True)

    return True


def midas_merge_genes(args):

    try:

        global global_args
        global_args = args

        global pool_of_samples
        global dict_of_species

        pool_of_samples = SamplePool(args.samples_list, args.midas_outdir, "genes")
        dict_of_species = pool_of_samples.select_species("genes", args)

        species_ids_of_interest = [sp.id for sp in dict_of_species.values()]
        assert species_ids_of_interest, f"No (specified) species pass the genome_coverage filter across samples, please adjust the genome_coverage or species_list"
        tsprint(species_ids_of_interest)

        pool_of_samples.create_dirs(["outdir", "tempdir"], args.debug)
        pool_of_samples.create_species_subdirs(species_ids_of_interest, "outdir", args.debug)
        pool_of_samples.create_species_subdirs(species_ids_of_interest, "tempdir", args.debug)
        pool_of_samples.write_summary_files(dict_of_species, "genes")

        # Download genes_info for every species in the restricted species profile.
        midas_db = MIDAS_DB(args.midas_db if args.midas_db else pool_of_samples.get_target_layout("midas_db_dir"))

        tsprint(f"CZ::design_chunks::start")
        arguments_list = design_chunks(midas_db, args.chunk_size)
        tsprint(f"CZ::design_chunks::finish")

        tsprint(f"CZ::fetch_clusters_map::start")
        for sp in dict_of_species.values():
            species_id = sp.id
            genes_info_file = midas_db.fetch_files("cluster_info", [species_id])[species_id]
            sp.fetch_clusters_map(genes_info_file, args.cluster_pid)
        tsprint(f"CZ::fetch_clusters_map::finish")


        # Merge copy_numbers, coverage and read counts across ALl the samples
        tsprint(f"CZ::multiprocessing_map::start")
        total_chunk_counts = sum((sp.num_of_genes_chunks for sp in dict_of_species.values()))
        num_cores = min(args.num_cores, total_chunk_counts)
        proc_flags = multiprocessing_map(process_one_chunk_of_genes, arguments_list, num_cores)
        assert all(s == "worked" for s in proc_flags)
        tsprint(f"CZ::multiprocessing_map::finish")

        tsprint(f"CZ::merge_all_chunks_per_species::start")
        for species_id in species_ids_of_interest:
            assert merge_all_chunks_per_species(species_id)
        tsprint(f"CZ::merge_all_chunks_per_species::start")

    except Exception as error:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error. Specify --debug flag to keep.")
            pool_of_samples.remove_dirs(["outdir", "tempdir"])
        raise error


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_merge_genes(args)
