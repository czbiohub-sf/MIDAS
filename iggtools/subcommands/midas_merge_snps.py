import json
import os
from collections import defaultdict
import multiprocessing

import numpy as np
from pysam import AlignmentFile  # pylint: disable=no-name-in-module
import Bio.SeqIO

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, num_physical_cores, command, InputStream, OutputStream, multithreading_hashmap, download_reference
from iggtools.params import outputs
from iggtools.models.uhgg import UHGG
from iggtools.common.samples import parse_species_profile, select_species
from iggtools.common.bowtie2 import build_bowtie2_db, bowtie2_align, samtools_index


DEFAULT_SAMPLE_COUNTS = 1
DEFAULT_SAMPLE_DEPTH = 5.0
DEFAULT_SAMPLE_COVERAGE = 0.4
DEFAULT_ALLELE_FREQ = 0.01
DEFAULT_SNP_TYPE = "bi"
DEFAULT_SITE_DEPTH = 1
DEFAULT_SITE_RATIO = 2.0
DEFAULT_SITE_PREV = 0.95
DEFAULT_SITE_MAF = 0.1

DECIMALS = ".3f"


snps_summary_schema = {
    "species_id": str,
    "genome_length": int,
    "covered_bases": int,
    "total_depth": int,
    "aligned_reads": int,
    "mapped_reads": int,
    "fraction_covered": float,
    "mean_coverage": float
}


def register_args(main_func):
    subparser = add_subcommand('midas_merge_snps', main_func, help='pooled-sample core genome SNPs calling')
    subparser.add_argument('outdir',
                           type=str,
                           help="""Path to directory to store results.  Subdirectory will be created for each species.""")
    subparser.add_argument('--sample_list',
                           dest='sample_list',
                           type=str,
                           required=True,
                           help=f"TSV file mapping sample name to midas_run_species.py output directories")

    # Species filters
    subparser.add_argument('--sample_counts',
                           dest='sample_counts', #min_samples
                           type=int,
                           metavar="INT",
                           default=DEFAULT_SAMPLE_COUNTS,
                           help=f"select species with >= MIN_SAMPLES ({DEFAULT_SAMPLE_COUNTS})")

    # Sample filters
    subparser.add_argument('--sample_depth',
                           dest='sample_depth',
                           type=float,
                           metavar="FLOAT",
                           default=DEFAULT_SAMPLE_DEPTH,
                           help=f"Minimum average read depth per sample ({DEFAULT_SAMPLE_DEPTH})")
    subparser.add_argument('--sample_coverage',
                           dest='sample_coverage', #fract_cov
                           type=float,
                           metavar="FLOAT",
                           default=DEFAULT_SAMPLE_COVERAGE,
                           help=f"Fraction of reference sites covered by at least 1 read ({DEFAULT_SAMPLE_COVERAGE})")

    # SNPs metadata
    subparser.add_argument('--allele_freq',
                           dest='allele_freq',
                           type=float,
                           metavar="FLOAT",
                           default=DEFAULT_ALLELE_FREQ,
                           help=f"Minimum frequency for calling an allele present ({DEFAULT_ALLELE_FREQ}), Values > 0.0 and < 0.5 are accepted.")
    subparser.add_argument('--snp_type',
                           type=str,
                           dest='snp_type',
                           default=DEFAULT_SNP_TYPE,
                           choices=['any','mono', 'bi', 'tri', 'quad'],
                           nargs='+',
                           help="""Specify one or more of the following:
                                    mono: keep sites with 1 allele > ALLELE_FREQ
                                    bi: keep sites with 2 alleles > ALLELE_FREQ (default)
                                    tri: keep sites with 3 alleles > ALLELE_FREQ
                                    quad: keep sites with 4 alleles > ALLELE_FREQ
                                    any: keep sites regardless of observed alleles""")

    # Site filters
    subparser.add_argument('--site_depth',
                           dest='site_depth',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_SITE_DEPTH,
                           help=f"Minimum number of reads mapped to genomic site ({DEFAULT_SITE_DEPTH})")
    subparser.add_argument('--site_ratio',
                           dest='site_ratio',
                           default=DEFAULT_SITE_RATIO,
                           type=float,
                           metavar="FLOAT",
                           help=f"Maximum ratio of site depth to genome depth ({DEFAULT_SITE_RATIO}).")
    subparser.add_argument('--site_prev',
                           dest='site_prev',
                           default=DEFAULT_SITE_PREV,
                           type=float,
                           metavar="FLOAT",
                           help=f"Minimum fraction of sample where genomic site satifying the site filters ({DEFAULT_SITE_PREV})")
    return main_func


def read_samples(sample_list):
    with InputStream(sample_list) as stream:
        samples = dict(select_from_tsv(stream, selected_columns=["sample_name", "midas_output_path"]))

    for sample_name in samples.keys():
        midas_output_path = samples[sample_name]
        assert os.path.exists(midas_output_path), f"MIDAS output directory {midas_output_path} for sample {sample_name} not exist."
        snps_summary = f"{midas_output_path}/snps/summary.txt"
        assert os.path.exists(species_profile), f"Missing MIDAS snps profile: {snps_summary}"
        samples[sample_name] = snps_summary
    return samples


def read_snps_summary(snps_summary_path):
    with InputStream(snps_summary_path) as instream:
        for record in select_from_tsv(instream, selected_columns=snps_summary_schema, result_structure=dict):
            yield record


def init_species(samples, args):
    "Store high quality sample-species pairs"

    species_samples = defaultdict(lambda: defaultdict(list))
    for sample_name, snps_summary_path in samples.items():
        # Read in the SNPs summary files
        for record in read_snps_summary(snps_summary_path):
            species_id = record["species_id"]
            mean_coverage = record["mean_coverage"]
            # Filter out low_coverage species-sample pair
            if mean_coverage < args.sample_depth:
                continue
            if  record["fraction_covered"] < args.sample_coverage:
                continue
            # Update the results
            if species_id in species_samples.keys():
                species_samples[species_id]["sample_counts"] += 1
            else:
                species_samples[species_id]["sample_counts"] = 1
            species_samples[species_id]["sample_name"].append(sample_name)
            species_samples[species_id]["sample_coverage"].append(mean_coverage)

    return species_samples


def sort_species(species_samples, sort_by="sample_counts"):
    # Order species by sample_counts in descentding order
    species_tuples = ((speices_id, species_samples[species_id][sort_by]) for species_id in species_samples.keys())
    sorted_species = sorted(species_tuples, key=lamdba x: x[1], reverse=True)
    return (_[0] for _ in sorted_species)


def filter_species(species_samples, args):
    "Pick subset of species to analyze"
    for species_id in sort_species(species_samples, "sample_counts"):
        if species_samples[species_id]["sample_counts"] < args.sample_counts:
            continue
    return species_samples


def select_species(samples, args):
    "Select species with a minimum number of high-coverage samples"
    species_samples = init_species(samples, args)
    species_samples = filter_species(species_samples, args)
    return species_samples


def midas_merge_snps(args):

    outdir = f"{args.outdir}/merged/snps"
    if args.debug and os.path.exists(outdir):
        tsprint(f"INFO:  Reusing existing output data in {outdir} according to --debug flag.")
    else:
        command(f"rm -rf {outdir}")
        command(f"mkdir -p {outdir}")

    # List samples and select species
    samples = read_samples(args.sample_list)
    species_samples = select_species(samples, args)

    # iggdb species_info => species_id	species_alt_id	species_name	representative_genome
    # iggdb genomes info => genome_id	genome_name	species_id
    # I want to leave create output directory not inside the select_species ...





@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_merge_snps(args)
