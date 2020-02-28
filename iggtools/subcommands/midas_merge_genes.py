import json
import os
from collections import defaultdict

from iggtools.models.pool import Pool, select_species

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, command, InputStream, OutputStream, select_from_tsv, multithreading_hashmap, multithreading_map, num_physical_cores, download_reference
from iggtools.params import outputs
from iggtools.params.schemas import genes_profile_schema, genes_info_schema, genes_schema, DECIMALS


DEFAULT_GENOME_DEPTH = 1.0
DEFAULT_SAMPLE_COUNTS = 1
DEFAULT_CLUSTER_ID = '95'
DEFAULT_MIN_COPY = 0.35


def register_args(main_func):
    subparser = add_subcommand('midas_merge_genes', main_func, help='metagenomic pan-genome profiling')

    subparser.add_argument('outdir',
                           type=str,
                           help="""Path to directory to store results.  Subdirectory will be created for each species.""")
    subparser.add_argument('--sample_list',
                           dest='sample_list',
                           type=str,
                           required=True,
                           help=f"TSV file mapping sample name to midas_run_species.py output directories")

    # Species and sample filters
    subparser.add_argument('--species_id',
                           dest='species_id',
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


def check_outdir(args):
    outdir = f"{args.outdir}/merged/genes"
    if args.debug and os.path.exists(outdir):
        tsprint(f"INFO:  Reusing existing output data in {outdir} according to --debug flag.")
    else:
        command(f"rm -rf {outdir}")
        command(f"mkdir -p {outdir}")

    paramstr = f"sd{args.genome_depth}.mc{args.min_copy}.cpid{args.cluster_pid}"
    tempdir = f"{outdir}/temp_{paramstr}"
    if args.debug and os.path.exists(tempdir):
        tsprint(f"INFO:  Reusing existing temp intermediate data in {tempdir} according to --debug flag.")
    else:
        command(f"rm -rf {tempdir}")
        command(f"mkdir -p {tempdir}")

    return (outdir, tempdir)


def pangenome_file(species_id, component):
    # s3://microbiome-igg/2.0/pangenomes/GUT_GENOMEDDDDDD/{genes.ffn, centroids.ffn, gene_info.txt}
    return f"{outputs.pangenomes}/{species_id}/{component}"


def read_cluster_map(gene_info_path, pid):
    genes_map = {}
    cols = ['gene_id', 'centroid_99', f"centroid_{pid}"]
    with InputStream(gene_info_path) as stream:
        for r in select_from_tsv(stream, selected_columns=cols, result_structure=dict):
            genes_map[r['centroid_99']] = r[f"centroid_{pid}"]
    return genes_map

def read_genes(path):
    with InputStream(path) as stream:
        for r in select_from_tsv(stream, selected_columns=genes_schema, result_structure=dict):
            yield r


def build_gene_matrices(sp, min_copy):
    """ Compute gene copy numbers for samples """
    for sample in sp.samples:
        sample.genes = {}
        for field, dtype in genes_info_schema.items():
            sample.genes[field] = defaultdict(dtype)

        path = f"{sample.dir}/genes/output/{sp.id}.genes.lz4"
        for r in read_genes(path):
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

    table_iterator = read_genes(midas_genes_dir)
    for r in table_iterator:
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


def per_species_worker(arguments):

    species_id, genes_info_path, species_samples, args, outdir = arguments
    sample_names = list(species_samples.keys())

    genes_map = read_cluster_map(genes_info_path, args.cluster_pid)

    accumulator = defaultdict(dict)
    for sample_index, sample_name in enumerate(sample_names):
        midas_genes_path = f"{species_samples[sample_name]}/genes/output/{species_id}.genes.lz4"
        assert midas_genes_path, f"Missing MIDAS genes output {midas_genes_path} for sample {sample_name}"

        ps_args = (genes_map, sample_index, midas_genes_path, len(species_samples))
        collect(accumulator, ps_args)

    for gene_id, copynum in accumulator["copynum"].items():
        accumulator["presabs"][gene_id] = [1 if cn > args.min_copy else 0 for cn in copynum]

    write(accumulator, sample_names, outdir, species_id)


def midas_merge_genes(args):

    outdir, tempdir = check_outdir(args)

    pool_of_samples = Pool(args.sample_list, "genes")
    list_of_species = select_species(pool_of_samples, args, "genes")

    for sp in list_of_species:
        print([sd.sample_name for sd in sp.samples])
        print(sp.samples_depth)
    exit(0)

    #list_of_species = select_species(args, "genes", outdir)

    # Download gene_info for every species in the restricted species profile.
    def download_genes_info(species_id):
        return download_reference(pangenome_file(species_id, "gene_info.txt.lz4"), f"{tempdir}/{species_id}")
    genes_info_files = multithreading_hashmap(download_genes_info, [sp.id for sp in list_of_species], num_threads=20)
    # Do we really need the hashmap ?? I don't think so

    # Collect copy_number, coverage and read counts across ALl the samples
    argument_list = []
    for sp in list_of_species:
        sp.write_summary("genes", outdir)
        species_id = sp.id
        species_samples = {sample.id: sample.dir for sample in sp.samples}
        argument_list.append((species_id, genes_info_files[species_id], species_samples, args, outdir))

    multithreading_map(per_species_worker, argument_list, num_threads=num_physical_cores)


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_merge_genes(args)
