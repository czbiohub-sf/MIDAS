import json
from collections import defaultdict

from iggtools.models.pool import Pool, select_species

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, command, InputStream, OutputStream, select_from_tsv, multiprocessing_map, num_physical_cores, download_reference
from iggtools.params import outputs
from iggtools.models.uhgg import UHGG
from iggtools.params.schemas import genes_summary_schema, genes_info_schema, genes_coverage_schema, format_data


DEFAULT_GENOME_DEPTH = 1.0
DEFAULT_SAMPLE_COUNTS = 1
DEFAULT_CLUSTER_ID = '95'
DEFAULT_MIN_COPY = 0.35

DEFAULT_CHUNK_SIZE = 5000


def read_cluster_map(packed_args):
    gene_info_path, pid = packed_args

    genes_map_dict = {}
    cols = ['gene_id', 'centroid_99', f"centroid_{pid}"]
    with InputStream(gene_info_path) as stream:
        for r in select_from_tsv(stream, selected_columns=cols, result_structure=dict):
            genes_map_dict[r['centroid_99']] = r[f"centroid_{pid}"]
    return genes_map_dict


def write_per_species(accumulator, species_id):

    global species_sliced_coverage_path
    global pool_of_samples

    sample_names = species_sliced_coverage_path[species_id]["species_samples"].keys()

    for file_type in genes_info_schema.keys():
        outfile = pool_of_samples.get_target_layout(f"genes_{file_type}", species_id)
        with OutputStream(outfile) as stream:
            stream.write("\t".join(["gene_id"] + sample_names) + "\n")
            for gene_id, gene_vals in accumulator[file_type].items():
                stream.write("\t".join(map(format_data, [gene_id] + gene_vals)) + "\n")


def collect(accumulator, my_args):

    species_id, sample_index, midas_genes_dir, total_sample_counts = my_args

    global species_sliced_coverage_path
    genes_map_dict = species_sliced_coverage_path[species_id]["genes_map"]

    with InputStream(midas_genes_dir) as stream:
        for r in select_from_tsv(stream, selected_columns=genes_info_schema, result_structure=dict):
            gene_id = genes_map_dict[r["gene_id"]]

            acc_copynum = accumulator["copynum"].get(gene_id)
            if not acc_copynum:
                acc_copynum = [0.0] * total_sample_counts
                accumulator["copynum"][gene_id] = acc_copynum
            acc_copynum[sample_index] += r["copy_number"]

            acc_depth = accumulator["depth"].get(gene_id)
            if not acc_depth:
                acc_depth = [0.0] * total_sample_counts
                accumulator["depth"][gene_id] = acc_depth
            acc_depth[sample_index] += r["coverage"]

            acc_reads = accumulator["reads"].get(gene_id)
            if not acc_reads:
                acc_reads = [0.0] * total_sample_counts
                accumulator["reads"][gene_id] = acc_reads
            acc_reads[sample_index] += r["count_reads"]


def per_species_worker(species_id):

    global pool_of_samples
    global global_args
    global species_sliced_coverage_path

    # Read in genes_info for each species
    #genes_info_path = pool_of_samples.get_target_layout("genes_info_file", species_id)
    #genes_map_dict = read_cluster_map(genes_info_path, global_args.cluster_pid)

    genes_map_dict = species_sliced_coverage_path[species_id]["genes_map"]
    species_samples = species_sliced_coverage_path[species_id]["species_samples"]
    sample_names = list(species_samples.keys())

    accumulator = defaultdict(dict)
    for sample_index, sample_name in enumerate(sample_names):
        midas_genes_path = species_samples[sample_name]
        my_args = (species_id, sample_index, midas_genes_path, genes_map_dict, len(species_samples))
        collect(accumulator, my_args)

    for gene_id, copynum in accumulator["copynum"].items():
        accumulator["presabs"][gene_id] = [1 if cn > global_args.min_copy else 0 for cn in copynum]

    write_per_species(accumulator, species_id)


def midas_merge_genes(args):

    global pool_of_samples

    pool_of_samples = Pool(args.samples_list, args.midas_outdir, "genes")
    list_of_species = select_species(pool_of_samples, "genes", args)

    # Create the output/temp (species) directory when we know what species to process
    species_ids_of_interest = [sp.id for sp in list_of_species]
    pool_of_samples.create_output_dir()
    pool_of_samples.create_species_subdir(species_ids_of_interest, "outdir", args.debug)
    pool_of_samples.create_species_subdir(species_ids_of_interest, "tempdir", args.debug)

    global global_args
    global_args = args

    global species_sliced_coverage_path
    species_sliced_coverage_path = defaultdict(dict)

    # Create species-to-process lookup table for each species
    local_toc = download_reference(outputs.genomes, pool_of_samples.get_target_layout("dbsdir"))
    db = UHGG(local_toc)
    # Download pan-genes-info for every species in the restricted species profile.
    genes_info_files = db.fetch_genes_info(species_ids_of_interest, pool_of_samples.get_target_layout("tempdir"))

    args_list = []
    for species_index, species_id in enumerate(species_ids_of_interest):
        #genes_info_path = pool_of_samples.get_target_layout("genes_info_file", species_id)
        genes_info_path = genes_info_files[species_index]
        args_list.append((genes_info_path, args.cluster_pid))
    list_of_species_genes_map = multiprocessing_map(read_cluster_map, args_list, num_physical_cores)

    for species_index, species_id in enumerate(species_ids_of_interest):
        species_sliced_coverage_path[species_id]["genes_map"] = list_of_species_genes_map[species_index]

    # Collect copy_number, coverage and read counts across ALl the samples
    argument_list = []
    for species in list_of_species:
        species_id = species.id

        # What is this? write_summary??
        species.write_summary("genes", pool_of_samples.get_target_layout("genes_summary", species_id))

        species_samples = dict()
        for sample in species.samples:
            sample_name = sample.sample_name
            midas_genes_path = sample.get_target_layout("genes_coverage", species_id)
            assert midas_genes_path, f"Missing MIDAS genes output {midas_genes_path} for sample {sample_name}"
            species_samples[sample.sample_name] = midas_genes_path
        species_sliced_coverage_path[species_id]["species_samples"] = species_samples

        argument_list.append(species_id)

    multiprocessing_map(per_species_worker, argument_list, num_physical_cores)


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


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_merge_genes(args)
