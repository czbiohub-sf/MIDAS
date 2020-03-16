import json
import os
from collections import defaultdict

from iggtools.models.pool import Pool, select_species

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, command, InputStream, OutputStream, select_from_tsv, multithreading_hashmap, multithreading_map, num_physical_cores, download_reference
from iggtools.params import outputs
from iggtools.params.schemas import genes_summary_schema, genes_info_schema, genes_coverage_schema, DECIMALS


DEFAULT_GENOME_DEPTH = 1.0
DEFAULT_SAMPLE_COUNTS = 1
DEFAULT_CLUSTER_ID = '95'
DEFAULT_MIN_COPY = 0.35

DEFAULT_CHUNK_SIZE = 5000


def read_cluster_map(gene_info_path, pid):
    genes_map = {}
    cols = ['gene_id', 'centroid_99', f"centroid_{pid}"]
    with InputStream(gene_info_path) as stream:
        for r in select_from_tsv(stream, selected_columns=cols, result_structure=dict):
            genes_map[r['centroid_99']] = r[f"centroid_{pid}"]
    return genes_map


def write(accumulator, sample_names, outdir, species_id):
    species_outdir = f"{outdir}/{species_id}"
    command(f"rm -rf {species_outdir}")
    command(f"mkdir -p {species_outdir}")

    for file_type in genes_info_schema.keys():
        outfile = f"{species_outdir}/genes_{file_type}.tsv.lz4"
        with OutputStream(outfile) as stream:
            stream.write("\t".join(["gene_id"] + sample_names) + "\n")
            for gene_id, gene_vals in accumulator[file_type].items():
                stream.write(f"{gene_id}\t" + "\t".join((str(format(val, DECIMALS)) for val in gene_vals)) + "\n")


def build_gene_matrices(sp, min_copy):
    """ Compute gene copy numbers for samples """
    for sample in sp.samples:
        sample.genes = {}
        for field, dtype in genes_info_schema.items():
            sample.genes[field] = defaultdict(dtype)

        path = f"{sample.dir}/genes/output/{sp.id}.genes.lz4"
        with InputStream(path) as stream:
            for r in select_from_tsv(stream, selected_columns=genes_schema, result_structure=dict):
                gene_id = sp.map[r["gene_id"]]
                sample.genes["copynum"][gene_id] += float(r['copy_number'])
                sample.genes["depth"][gene_id] += float(r["coverage"])
                sample.genes["reads"][gene_id] += int(r["count_reads"])

    for sample in sp.samples:
        for gene_id, copynum in sample.genes["copynum"].items():
            if copynum >= min_copy:
                sample.genes["presabs"][gene_id] = 1
            else:
                sample.genes["presabs"][gene_id] = 0


def write_gene_matrices(sp, outdir):
    """ Write pangenome matrices to file """

    sample_names = [s.id for s in sp.samples]
    # This doesn't make sense to me ...
    genes = sorted(sp.samples[0].genes["depth"])

    for file_type in genes_info_schema.keys():
        outfile = f"{outdir}/genes_{file_type}.tsv.lz4"
        with OutputStream(outfile) as stream:
            stream.write("\t".join(["gene_id"] + sample_names) + "\n")
            record = []
            for gene_id in genes:
                for sample in sp.samples:
                    record.append(sample.genes[file_type][gene_id])
            stream.write("\t".join(map(str, [gene_id] + record)) + "\n")


def collect(accumulator, ps_args):

    genes_map, sample_index, midas_genes_dir, total_sample_counts = ps_args

    with InputStream(midas_genes_dir) as stream:
        for r in select_from_tsv(stream, selected_columns=genes_schema, result_structure=dict):
            gene_id = genes_map[r["gene_id"]]

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


def per_species_worker(arguments):
#def compute_chunk_gene_coverage(arguments):

    global pool_of_samples
    global global_args

    #species_id, genes_info_path, species_samples, args, outdir = arguments
    species_id, species_samples, args, outdir = arguments
    sample_names = list(species_samples.keys())

    # Read in genes_info for each species
    genes_info_path = pool_of_samples.get_target_layout("genes_info_file", species_id)
    genes_map_dict = read_cluster_map(genes_info_path, global_args.cluster_pid)

    accumulator = defaultdict(dict)
    for sample_index, sample_name in enumerate(sample_names):
        #midas_genes_path = f"{species_samples[sample_name]}/genes/output/{species_id}.genes.lz4"
        midas_genes_path = sample.get_target_layout("genes_coverage", species_id)
        assert midas_genes_path, f"Missing MIDAS genes output {midas_genes_path} for sample {sample_name}"

        ps_args = (genes_map_dict, sample_index, midas_genes_path, len(species_samples))
        collect(accumulator, ps_args)

    for gene_id, copynum in accumulator["copynum"].items():
        accumulator["presabs"][gene_id] = [1 if cn > args.min_copy else 0 for cn in copynum]

    write(accumulator, sample_names, outdir, species_id)


def process_chunk_of_genes(packed_args):

    global semaphore_for_species
    global list_of_species
    global species_sliced_coverage_path

    if packed_args[1] == -1:
        species_id = packed_args[0]
        number_of_chunks = len(species_sliced_coverage_path[species_id])

    species_id, chunk_id, contig_id, contig_start, contig_end, samples_depth, samples_snps_pileup, total_samples_count = packed_args
    samples_depth
    samples_snps_pileup
    total_samples_count
    try:
        accumulator = defaultdict(dict)
        for sample_index in range(total_samples_count):
            # how to find Sample object by sample_index
            midas_genes_path = sample.get_target_layout("genes_coverage", species_id)

            ps_args = (genes_map_dict, sample_index, midas_genes_path, len(species_samples))
            collect(accumulator, ps_args)


    finally:
        semaphore_for_species[species_id].release() # no deadlock



def midas_merge_genes(args):

    global pool_of_samples
    global list_of_species

    pool_of_samples = Pool(args.samples_list, args.midas_outdir, "genes")
    list_of_species = select_species(pool_of_samples, "genes", args)

    # Create the output/temp (species) directory when we know what species to process
    species_ids_of_interest = [sp.id for sp in list_of_species]
    pool_of_samples.create_output_dir()
    pool_of_samples.create_species_subdir(species_ids_of_interest, "outdir", args.debug)
    pool_of_samples.create_species_subdir(species_ids_of_interest, "tempdir", args.debug)


    global global_args
    global_args = args

    # Create species-to-process lookup table for each species
    local_toc = download_reference(outputs.genomes, pool_of_samples.get_target_layout("dbsdir"))
    db = UHGG(local_toc)
    # Download pan-genes-info for every species in the restricted species profile.
    genes_info_files = db.fetch_genes_info(species_ids_of_interest, pool_of_samples.get_target_layout("tempdir"))

    # Collect copy_number, coverage and read counts across ALl the samples

    global species_sliced_coverage_path
    species_sliced_coverage_path = defaultdict(dict)

    argument_list = []
    for species in list_of_species:
        species.write_summary("genes", outdir)
        species_id = species.id
        species_sliced_coverage_path[species_id]["sample"] =

        for sample in species.samples:
            midas_genes_path = sample.get_target_layout("genes_coverage", species_id)
            sample.sample_name

        species_samples = {sample.id: sample.dir for sample in species.samples}

        for sample_index in range(total_samples_count):
            # how to find Sample object by sample_index



        argument_list.append((species_id, species_samples, args, outdir))

    multithreading_map(per_species_worker, argument_list, num_threads=num_physical_cores)


def register_args(main_func):
    subparser = add_subcommand('midas_merge_genes', main_func, help='metagenomic pan-genome profiling')

    subparser.add_argument('midas_outdir',
                           type=str,
                           help="""Path to directory to store results.  Subdirectory will be created for each species.""")
    subparser.add_argument('--sample_list',
                           dest='sample_list',
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
