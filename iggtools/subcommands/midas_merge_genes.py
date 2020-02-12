import json
import os
from collections import defaultdict
import numpy as np
import Bio.SeqIO
from pysam import AlignmentFile  # pylint: disable=no-name-in-module

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, command, InputStream, OutputStream, select_from_tsv, multithreading_hashmap, download_reference
from iggtools.params import outputs
from iggtools.common.samples import parse_species_profile, select_species
from iggtools.common.bowtie2 import build_bowtie2_db, bowtie2_align

CLUSTERING_PERCENTS = [99, 95, 90, 85, 80, 75]

DEFAULT_ALN_COV = 0.75
DEFAULT_SPECIES_COVERAGE = 3.0
DEFAULT_ALN_MAPID = 94.0
DEFAULT_ALN_READQ = 20
DEFAULT_ALN_MAPQ = 0

DECIMALS = ".6f"


genes_summary_schema = {
    "species_id": str,
    "pangenome_size": int,
    "covered_genes": int,
    "fraction_covered": float,
    "mean_coverage": float,
    "aligned_reads": int,
    "mapped_reads": int
}


def register_args(main_func):
    subparser = add_subcommand('midas_run_genes', main_func, help='metagenomic pan-genome profiling')
    subparser.add_argument('outdir',
                           type=str,
                           help="""Path to directory to store results.  Name should correspond to unique sample identifier.""")
    subparser.add_argument('-1',
                           dest='r1',
                           required=True,
                           help="FASTA/FASTQ file containing 1st mate if using paired-end reads.  Otherwise FASTA/FASTQ containing unpaired reads.")
    subparser.add_argument('-2',
                           dest='r2',
                           help="FASTA/FASTQ file containing 2nd mate if using paired-end reads.")
    subparser.add_argument('--max_reads',
                           dest='max_reads',
                           type=int,
                           metavar="INT",
                           help=f"Number of reads to use from input file(s).  (All)")
    subparser.add_argument('--species_cov',
                           type=float,
                           dest='species_cov',
                           metavar='FLOAT',
                           default=DEFAULT_SPECIES_COVERAGE,
                           help=f"Include species with >X coverage ({DEFAULT_SPECIES_COVERAGE})")
    if False:
        # This is unused.
        subparser.add_argument('--species_topn',
                               type=int,
                               dest='species_topn',
                               metavar='INT',
                               help='Include top N most abundant species')

    #  Alignment flags (bowtie, or postprocessing)
    subparser.add_argument('--aln_cov',
                           dest='aln_cov',
                           default=DEFAULT_ALN_COV,
                           type=float,
                           metavar="FLOAT",
                           help=f"Discard reads with alignment coverage < ALN_COV ({DEFAULT_ALN_COV}).  Values between 0-1 accepted.")
    subparser.add_argument('--aln_readq',
                           dest='aln_readq',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_ALN_READQ,
                           help=f"Discard reads with mean quality < READQ ({DEFAULT_ALN_READQ})")
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
                           help='Global/local read alignment (local, corresponds to the bowtie2 --local).  Global corresponds to the bowtie2 default --end-to-end.')
    subparser.add_argument('--aln_interleaved',
                           action='store_true',
                           default=False,
                           help='FASTA/FASTQ file in -1 are paired and contain forward AND reverse reads')

    return main_func


class Sample:
    """
    base class for sample, to keep track of the midas species/snps/genes path
    """
    def __init__(self, midas_outdir, sample_name, dbtype):
        self.dir = midas_outdir
        self.id = sample_name
        self.info = self.read_info(dbtype)

    def read_info(self, dbtype):

        assert os.path.exists(self.dir), f"Invalid MIDAS output {self.dir} for {self.id}"

        summary_path = f"{self.dir}/{dbtype}/summary.txt"
        assert os.path.exists(summary_path), f"Missing MIDAS {dbtype} summary for {self.id}"

        if dbtype == "snps":
            schema = snps_summary_schema
        if dbtype = "genes":
            schema = genes_summary_schema

        info = defaultdict()
        with InputStream(summary_path) as stream:
            for row in select_from_tsv(stream, selected_columns=schema, result_structure=dict):
                info[row["species_id"]] = row
        return info


class Species:
    """
    Base class for species
    """
    def __init__(self, id, species_info, genome_info):
        self.species_id = id
        self.samples = []

    def fetch_sample_depth(self):
        self.sample_depth = [sample.info[self.species_id]["mean_coverage"] for sample in self.samples]


def read_samples(sample_list, dbtype):
    """
    Input: read in Table-of-Content: sample_name\tpath/to/midas/snps_dir

    Output: return dict of samples: local file path of midas_outdir
    """
    samples = []
    with InputStream(sample_list) as stream:
        for row in select_from_tsv(stream, selected_columns=["sample_name", "midas_output_path"]):
            sample = Sample(row["midas_output_path"], row["sample_name"], dbtype)
            samples.append(sample)
    return samples


def filter_sample_species(info, species_id, args, dbtype):
    """
    Select high quality sample-species pairs based on the genome_coverage and
    genome_depth
    """

    if info['mean_coverage'] < args.sample_depth:
        return False # skip low-coverage <species, sample>
    if (args.species_id and species_id not in args.species_id.split(",")):
        return False # skip unspeficied species
    if (dbtype == "snps" and info['fraction_covered'] < args.genome_coverage):
        return False # skip low prevalent <species, sample>
    return True


def init_species(samples, args, dbtype):
    species = {}
    for sample in samples:
        for species_info in sample.info.values():
            species_id = species_info["species_id"]
            if species_id not in species:
                species[species_id] = Species(species_id)

            if filter_samples_species(species_info, species_id, args, dbtype):
                species[species_id].samples.append(sample)
    return list(species.values())


def sort_species(species):
    """
    Sort species by sample_counts in descending order
    Output: list of sorted Species objects
    """

    species_sorted = sorted(((sp, len(sp.samples)) for sp in species), key=lambda x: x[1], reverse = True)

    return [sp[0] for sp in species_sorted]


def filter_species(species, args):
    """
    Filter out low prevalent species using sample_counts cutoff.
    """

    keep = []
    species_genes_summary = []

    for sp in sorted_species(species):
        # Pick subset of species to analyze
        sp.sample_counts = len(sp.samples)

        if sp.sample_counts < args.min_samples:
            continue # skip low prevalent species
        if (args.max_samples and sp.sample_counts >= args.max_samples):
            continue # skip species with too many samples
        if (args.max_species and len(keep) > args.max_species):
            continue

        sp.fetch_sample_depth()
        keep.append(sp)

        for sample in sp.samples:
            species_genes_summary.append([sp.species_id, sample.id] + list(sample.info[sp.species_id].values())[1:])

    return (keep, species_genes_summary)


def write_summary(species_genes_summary, outdir):
    """ Write summary file for samples """
    genes_summary_dir = f"{outdir}/genes_summary.tsv.lz4"
    with OutputStream(genes_summary_dir) as stream:
        stream.write("\t".join(["species_id", "sample_id"] + list(genes_summary_schema.keys())[1:]) + "\n")
        for record in species_genes_summary:
            stream.write("\t".join(map(str, record)) + "\n")


def select_species(args, dbtype):
    # Read in list of samples, return [Sample1, Sample2, .., Samplen]
    samples = init_samples(args.sample_list, dbtype)
    species = init_species(samples, args, dbtype)
    species, genes_summary = filter_species(species, args)
    write_summary(genes_summary, outdir)

    return species


def pangenome_file(species_id, component):
    # s3://microbiome-igg/2.0/pangenomes/GUT_GENOMEDDDDDD/{genes.ffn, centroids.ffn, gene_info.txt}
    return f"{outputs.pangenomes}/{species_id}/{component}"


def read_cluster_map(sp, pid, tempdir):
    sp.map = {}
    header = ['gene_id'] + [f"centroid_{pid}" for pid in CLUSTERING_PERCENTS]

    with InputStream(sp.gene_info) as stream:
        for row in select_from_tsv(stream, selected_columns=header, result_structure=dict):
            sp.map[r['centroid_99']] = r[f"centroid_{pid}}"]


def midas_merge_genes(args):

    outdir = f"{args.outdir}/merged/genes"
    if args.debug and os.path.exists(outdir):
        tsprint(f"INFO:  Reusing existing output data in {outdir} according to --debug flag.")
    else:
        command(f"rm -rf {outdir}")
        command(f"mkdir -p {outdir}")

    paramstr = f"sd{args.site_depth}"
    tempdir = f"{outdir}/temp_{paramstr}"
    if args.debug and os.path.exists(tempdir):
        tsprint(f"INFO:  Reusing existing temp intermediate data in {tempdir} according to --debug flag.")
    else:
        command(f"rm -rf {tempdir}")
        command(f"mkdir -p {tempdir}")


    species_list = select_species(args, dbtype = "genes", outdir)

    def download_genes_info(species_id):
        return download_reference(pangenome_file(species_id, "gene_info.txt.lz4"), f"{tempdir}/{species_id}")

    # Download gene_info for every species in the restricted species profile.
    genes_info_files = multithreading_hashmap(download_genes_info, [sp.speices_id for sp in species_list], num_threads=20)
    for sp in species_list:
        sp.gene_info = genes_info_file[sp.species_id]


    # Accumulate read_counts and sample_counts across ALL the sample passing the genome filters;
    # and at the same time remember <site, sample>'s A, C, G, T read counts.
    argument_list = []
    for species_id in list(species_samples.keys()):
        argument_list.append((species_id, species_samples[species_id], species_contigs[species_id], samples_midas, args, outdir))

    multithreading_map(per_species_worker, argument_list, num_threads=num_physical_cores)



@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_run_genes(args)
