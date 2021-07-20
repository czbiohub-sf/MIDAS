#!/usr/bin/env python3
import os
import json
from collections import defaultdict
import multiprocessing
from itertools import repeat

from iggtools.models.samplepool import SamplePool
from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, InputStream, OutputStream, select_from_tsv
from iggtools.models.midasdb import MIDAS_DB
from iggtools.params.schemas import genes_info_schema, genes_coverage_schema, format_data, fetch_default_genome_depth, DECIMALS6


DEFAULT_GENOME_DEPTH = fetch_default_genome_depth("genes")
DEFAULT_SAMPLE_COUNTS = 1
DEFAULT_CLUSTER_ID = '95'
DEFAULT_MIN_COPY = 0.35


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


def collect(accumulator, my_args):

    species_id, sample_index, genes_coverage_path, pid, total_samples_count, cluster_info = my_args

    with InputStream(genes_coverage_path) as stream:
        for r in select_from_tsv(stream, selected_columns=genes_coverage_schema, result_structure=dict):

            gene_id = cluster_info[r["gene_id"]][f"centroid_{pid}"]

            acc_copynum = accumulator["copynum"].get(gene_id)
            if not acc_copynum:
                acc_copynum = [0.0] * total_samples_count
                accumulator["copynum"][gene_id] = acc_copynum
            acc_copynum[sample_index] += r["copy_number"]

            acc_depth = accumulator["depth"].get(gene_id)
            if not acc_depth:
                acc_depth = [0.0] * total_samples_count
                accumulator["depth"][gene_id] = acc_depth
            acc_depth[sample_index] += r["mean_coverage"]

            acc_reads = accumulator["reads"].get(gene_id)
            if not acc_reads:
                acc_reads = [0.0] * total_samples_count
                accumulator["reads"][gene_id] = acc_reads
            acc_reads[sample_index] += r["mapped_reads"]


def build_gene_matrices(sp, pid, args_mincopy):

    species_id = sp.id
    total_samples_count = sp.samples_count
    cluster_info = sp.cluster_info

    tsprint(f"    CZ::process_one_chunk_of_genes::{species_id}::start collect_sample_by_sample")

    accumulator = defaultdict(dict)
    for sample_index, sample in enumerate(sp.list_of_samples):
        tsprint(f"    CZ::process_one_chunk_of_genes::{species_id}-{sample_index}::start collect_one_sample")
        genes_coverage_path = sample.get_target_layout("genes_coverage", species_id)

        my_args = (species_id, sample_index, genes_coverage_path, pid, total_samples_count, cluster_info)
        collect(accumulator, my_args)
        tsprint(f"    CZ::process_one_chunk_of_genes::{species_id}-{sample_index}::finish collect_one_sample")

    tsprint(f"    CZ::process_one_chunk_of_genes::{species_id}::finish collect_sample_by_sample")

    # Second pass: infer presence absence based on copy number
    for gene_id, copynum in accumulator["copynum"].items():
        accumulator["presabs"][gene_id] = [1 if cn >= args_mincopy else 0 for cn in copynum]

    return accumulator


def write_gene_matrices(accumulator, pool_of_samples, species_id, samples_names):
    for file_type in list(genes_info_schema.keys()):
        outfile = pool_of_samples.get_target_layout(f"genes_{file_type}", species_id)
        with OutputStream(outfile) as stream:
            stream.write("gene_id\t" + "\t".join(samples_names) + "\n")
            for gene_id, gene_vals in accumulator[file_type].items():
                stream.write(f"{gene_id}\t" + "\t".join(map(format_data, gene_vals, repeat(DECIMALS6, len(gene_vals)))) + "\n")
    return True


def midas_merge_genes(args):

    try:
        pool_of_samples = SamplePool(args.samples_list, args.midas_outdir, "genes")
        dict_of_species = pool_of_samples.select_species("genes", args)

        species_ids_of_interest = [sp.id for sp in dict_of_species.values()]
        species_counts = len(species_ids_of_interest)
        assert species_ids_of_interest, f"No (specified) species pass the genome_coverage filter across samples, please adjust the genome_coverage or species_list"
        tsprint(species_ids_of_interest)

        pool_of_samples.create_dirs(["outdir"], args.debug)
        pool_of_samples.create_species_subdirs(species_ids_of_interest, "outdir", args.debug)
        pool_of_samples.write_summary_files(dict_of_species, "genes")

        # Download genes_info for every species in the restricted species profile.
        midas_db = MIDAS_DB(args.midas_db if args.midas_db else pool_of_samples.get_target_layout("midas_db_dir"), 1)

        # Merge copy_numbers, coverage and read counts across ALl the samples
        tsprint(f"CZ::build_and_write_gene_matrices::start")

        for species_id in species_ids_of_interest:
            sp = dict_of_species[species_id]
            sp.get_cluster_info(midas_db)

            tsprint(f"  CZ::midas_merge_genes::{species_id}::start build_gene_matrices")
            accumulator = build_gene_matrices(sp, args.cluster_pid, args.min_copy)
            tsprint(f"  CZ::midas_merge_genes::{species_id}::finish build_gene_matrices")

            tsprint(f"  CZ::midas_merge_genes::{species_id}::start write_gene_matrices")
            samples_names = sp.fetch_samples_names()
            assert write_gene_matrices(accumulator, pool_of_samples, species_id, samples_names)
            tsprint(f"  CZ::midas_merge_genes::{species_id}::finish write_gene_matrices")

        tsprint(f"CZ::build_and_write_gene_matrices::finish")

    except Exception as error:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error. Specify --debug flag to keep.")
            pool_of_samples.remove_dirs(["outdir", "tempdir"])
        raise error


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_merge_genes(args)
