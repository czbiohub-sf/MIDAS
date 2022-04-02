#!/usr/bin/env python3
import os
import json
from collections import defaultdict
from itertools import repeat
from math import ceil

from midas2.models.samplepool import SamplePool
from midas2.common.argparser import add_subcommand
from midas2.common.utils import tsprint, InputStream, OutputStream, select_from_tsv, multiprocessing_map
from midas2.models.midasdb import MIDAS_DB
from midas2.params.schemas import genes_info_schema, genes_coverage_schema, format_data, DECIMALS6
from midas2.models.species import scan_cluster_info
from midas2.params.inputs import MIDASDB_NAMES


DEFAULT_GENOME_DEPTH = 1.0
DEFAULT_SAMPLE_COUNTS = 1
DEFAULT_CLUSTER_ID = '95'
DEFAULT_MIN_COPY = 0.35
DEFAULT_NUM_CORES = 4


def register_args(main_func):
    subparser = add_subcommand('merge_genes', main_func, help='metagenomic pan-genome profiling')

    subparser.add_argument('midas_outdir',
                           type=str,
                           help="""Path to directory to store results.  Subdirectory will be created for each species.""")
    subparser.add_argument('--samples_list',
                           dest='samples_list',
                           type=str,
                           required=True,
                           help=f"TSV file mapping sample name to run_species.py output directories")

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
    subparser.add_argument('--midas_db',
                           dest='midas_db',
                           type=str,
                           metavar="CHAR",
                           help=f"local MIDAS DB which mirrors the s3 IGG db")
    subparser.add_argument('--num_cores',
                           dest='num_cores',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_NUM_CORES,
                           help=f"Number of physical cores to use ({DEFAULT_NUM_CORES})")

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


def process(list_of_species):
    for species_id in list_of_species:
        accumulator = build_gene_matrices(species_id)
        assert write_gene_matrices(accumulator, species_id)
    return "worked"


def build_gene_matrices(species_id):
    global global_args
    global dict_of_species

    min_copy = global_args.min_copy
    sp = dict_of_species[species_id]

    # First pass: accumulate the gene matrix sample by sample
    accumulator = defaultdict(dict)
    for sample_index, sample in enumerate(sp.list_of_samples):
        genes_coverage_fp = sample.get_target_layout("genes_coverage", species_id)
        my_args = (species_id, sample_index, genes_coverage_fp)
        collect(accumulator, my_args)

    # Second pass: infer presence absence based on copy number
    for gene_id, copynum in accumulator["copynum"].items():
        accumulator["presabs"][gene_id] = [1 if cn >= min_copy else 0 for cn in copynum]

    return accumulator


def collect(accumulator, my_args):
    # Merge copy_numbers, coverage and read counts across ALl the samples

    global global_args
    global pool_of_samples
    global dict_of_species

    pid = global_args.cluster_pid

    species_id, sample_index, genes_coverage_fp = my_args

    sp = dict_of_species[species_id]
    total_samples_count = sp.samples_count
    cluster_info = scan_cluster_info(sp.cluster_info_fp)

    with InputStream(genes_coverage_fp) as stream:
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


def write_gene_matrices(accumulator, species_id):
    global dict_of_species
    global pool_of_samples

    samples_names = dict_of_species[species_id].fetch_samples_names()

    for file_type in list(genes_info_schema.keys()):
        outfile = pool_of_samples.get_target_layout(f"genes_{file_type}", species_id)
        with OutputStream(outfile) as stream:
            stream.write("gene_id\t" + "\t".join(samples_names) + "\n")
            for gene_id, gene_vals in accumulator[file_type].items():
                stream.write(f"{gene_id}\t" + "\t".join(map(format_data, gene_vals, repeat(DECIMALS6, len(gene_vals)))) + "\n")
    return True


def merge_genes(args):

    try:
        global global_args
        global_args = args

        global pool_of_samples
        global dict_of_species

        pool_of_samples = SamplePool(args.samples_list, args.midas_outdir, "genes")
        dict_of_species = pool_of_samples.select_species("genes", args)

        species_ids_of_interest = [sp.id for sp in dict_of_species.values()]
        species_counts = len(species_ids_of_interest)

        assert species_ids_of_interest, f"No (specified) species pass the genome_coverage filter across samples, please adjust the genome_coverage or species_list"
        tsprint(f"{species_counts} species pass the filter")

        pool_of_samples.create_dirs(["outdir"], args.debug)
        pool_of_samples.create_species_subdirs(species_ids_of_interest, "outdir", args.debug, quiet=True)
        pool_of_samples.write_summary_files(dict_of_species, "genes")

        # Download genes_info for every species in the restricted species profile.
        def chunkify(L, n):
            return [L[x: x+n] for x in range(0, len(L), n)]
        chunk_size = ceil(species_counts / args.num_cores)

        midas_db = MIDAS_DB(os.path.abspath(args.midasdb_dir), args.midasdb_name, args.num_cores)
        midas_db.fetch_files("pangenome_cluster_info", species_ids_of_interest)

        for species_id in species_ids_of_interest:
            sp = dict_of_species[species_id]
            sp.get_cluster_info_fp(midas_db)

        chunk_lists = chunkify(species_ids_of_interest, chunk_size)
        proc_flags = multiprocessing_map(process, chunk_lists, args.num_cores)
        assert all(s == "worked" for s in proc_flags)

    except Exception as error:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error. Specify --debug flag to keep.")
            pool_of_samples.remove_dirs(["outdir"])
        raise error


@register_args
def main(args):
    tsprint(f"Merge pangenome coverage in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    merge_genes(args)
