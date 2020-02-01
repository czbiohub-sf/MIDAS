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
DEFAULT_GENOME_DEPTH = 5.0
DEFAULT_GENOME_COVERAGE = 0.4
DEFAULT_ALLELE_FREQ = 0.01
DEFAULT_SNP_TYPE = "bi"
DEFAULT_SITE_DEPTH = 1
DEFAULT_SITE_RATIO = 2.0
DEFAULT_SITE_PREV = 0.95
DEFAULT_SITE_MAF = 0.1

DEBUG_MAX_LINES = 1000 * 1000

DECIMALS = ".3f"


## TODO: I think this is genome_stats, where is my contigs_stats?
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


snps_pileup_schema = {
    "ref_id": str,
    "ref_pos": int,
    "ref_allele": str,
    "depth": int,
    "count_a": int,
    "count_c": int,
    "count_g": int,
    "count_t": int
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
    subparser.add_argument('--genome_depth',
                           dest='genome_depth',
                           type=float,
                           metavar="FLOAT",
                           default=DEFAULT_genome_depth,
                           help=f"Minimum average read depth per sample ({DEFAULT_genome_depth})")
    subparser.add_argument('--genome_coverage',
                           dest='genome_coverage', #fract_cov
                           type=float,
                           metavar="FLOAT",
                           default=DEFAULT_genome_coverage,
                           help=f"Fraction of reference sites covered by at least 1 read ({DEFAULT_genome_coverage})")

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
    "Return dict of sample: midas_snps_summary_path"
    with InputStream(sample_list) as stream:
        samples = dict(select_from_tsv(stream, selected_columns=["sample_names", "midas_output_path"]))

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


def sort_species(species_samples, sort_by="sample_counts"):
    "Order species by sample_counts in descentding order"
    species_tuples = ((speices_id, species_samples[species_id][sort_by]) for species_id in species_samples.keys())
    sorted_species = sorted(species_tuples, key=lamdba x: x[1], reverse=True)
    return (_[0] for _ in sorted_species)


def select_species(samples, args):
    """
    Select high quality sample-species pairs
    For the given list of samples, read in snps_summary and filter species meeting
    the vertical coverage (genome_depth) and horizontal coverage (genome_coverage).
    For each species, we a dict with two lists: sample_name and genome_coverage.

    Further filter out species with few sample_counts.
    """
    species_samples = defaultdict(lambda: defaultdict(list))
    for sample_name, snps_summary_path in samples.items():
        # Read in the SNPs summary files
        for record in read_snps_summary(snps_summary_path):
            species_id = record["species_id"]
            mean_coverage = record["mean_coverage"]
            # Filter out low_coverage species-sample pair
            if mean_coverage < args.genome_depth:
                continue
            if record["fraction_covered"] < args.genome_coverage:
                continue
            # Update the results
            if species_id in species_samples.keys():
                species_samples[species_id]["sample_counts"] += 1
            else:
                species_samples[species_id]["sample_counts"] = 1

            species_samples[species_id]["sample_names"].append(sample_name)
            species_samples[species_id]["genome_coverage"].append(mean_coverage)

    for species_id in sort_species(species_samples, "sample_counts"):
        if species_samples[species_id]["sample_counts"] < args.sample_counts:
            continue

    return species_samples


def tsv_rows_slice(path, thread_nums, thread_id):
    with InputStream(path) as stream:
        for line in stream:
            yield next(stream).rstrip("\n").split("\t")
            for line in stream:
                if int(line[4:6] % thread_nums == thread_id):
                    yield line.rstrip("\n").split("\t")


def read_snps_pileup(snps_pile_path):
    with InputStream(snps_pile_path) as instream:
        for record in select_from_tsv(instream, selected_columns=snps_pileup_schema, result_structure=dict):
            yield record


def accumulate(accumulator, ps_args):
    """
    Process one sample
    """
    species_id, sample_index, sample_name, sample_counts, midas_snps_dir, genome_coverage = ps_args

    table_iterator = read_snps_pileup(midas_snps_dir)

    # Output column indices
    s_A, s_C, s_G, s_T, s_sample_count, s_scA, s_scC, s_scG, s_scT = range(9)

    for line, row in enumerate(table_iterator):
        if line == DEBUG_MAX_LINES:
            break

        # Unpack frequently accessed columns
        contig_id = row["ref_id"]
        ref_pos = row["ref_pos"]
        ref_allele = row["ref_allele"]
        depth = row["depth"]

        # Computae derived columns
        site_id = f"{contig_id}|{ref_pos}|{ref_allele}"

        A, C, G, T = row["count_a"], row["count_c"], row["count_g"], row["count_t"]
        number_alleles = 0
        nonzero_allele_index = 4
        sc_ACGT = [0, 0, 0, 0] # sample count for A C G T
        for i, nt_count in enumerate((A, C, G, T)):
            if nt_count / depth >= args.allele_freq:
                number_alleles += 1
                nonzero_allele_index = i
                sc_ACGT[i] = 1

        nz_allele = "ACGTN"[nonzero_allele_index]
        nz_allele_count = (A, C, G, T)[nonzero_allele_index]
        nz_allele_freq = nz_allele_count / depth
        site_ratio = depth / genome_coverage #<- this used to be contigs_coverage

        # Per-sample site filters
        if depth < args.site_depth:
            continue
        if site_ratio > args.site_ratio:
            continue
        if number_alleles > 2: #<- TODO: refer to MIDAS to handle more than 2 SNPs
            continue

        # TODO: for each site, should I really put a sample_names_lists (even of sample index)?

        # Aggragate
        species_acc = accumulator[species_id]
        acc = species_acc.get(site_id)
        if acc:
            acc[s_A] += A
            acc[s_C] += C
            acc[s_G] += G
            acc[s_T] += T
            acc[s_sample_count] += 1
            acc[s_scA] += sc_ACGT[0]
            acc[s_scC] += sc_ACGT[1]
            acc[s_scG] += sc_ACGT[2]
            acc[s_scT] += sc_ACGT[3]
        else:
            acc = [A, C, G, T, 1, sc_ACGT[0], sc_ACGT[1], sc_ACGT[2], sc_ACGT[3]] + ([('N', 0)] * sample_counts)
            species_acc[site_id] = acc

        # Slice by samples
        # This isn't being accumulated across samples;
        # we are just remembering the value from each sample.
        assert acc[9 + sample_index] == ('N', 0) and nz_allele != 'N'
        acc[9 + sample_index] = (nz_allele, nz_allele_freq)
        ## TODO: MIDAS seems to keep record of the minor_allele_frequency.


def write(accumulator, sample_names, outdir):
    "Write the matrix together with desired matrix"

    for species_id, species_acc in accumulator.items():

        output_sites = f"{outdir}/spid.{species_id}/snps_freq.tsv"

        with OutputStream(output_sites) as out_stream:
            out_sites.write("site_id\tA\tC\tG\tT\tsample_count\tscA\tscC\tscG\tscT\t")
            out_sites.write("\t".join(["major_allele", "minor_allele"] + sample_brief_names) + "\n")

            for site_id, site_info in genome_acc.items():
                A, C, G, T, sample_count, scA, scC, scG, scT = site_info[:9]
                depth = A + C + G + T
                all_alleles = ((A, 'A'), (C, 'C'), (G, 'G'), (T, 'T'))
                alleles_above_cutoff = tuple(al for al in all_alleles if al[0] / depth >= param.MIN_ALLELE_FREQUECY_ACROSS_SAMPLES)

                # Keep only bi-allelic and mono-allelic sites.
                if 1 <= len(alleles_above_cutoff) <= 2:
                    # In the event of a tie -- biallelic site with 50/50 freq split -- the allele declared major is
                    # the one that comes later in the "ACGT" lexicographic order.
                    alleles_above_cutoff = sorted(alleles_above_cutoff, reverse=True)
                    major_allele = alleles_above_cutoff[0][1]
                    minor_allele = alleles_above_cutoff[-1][1]  # for mono-allelic sites, same as major allele
                    out_sites.write(f"{site_id}\t{A}\t{C}\t{G}\t{T}\t{sample_count}\t{scA}\t{scC}\t{scG}\t{scT}\t")
                    ## I vaguely remembered Katie also suggested doing minor_allele_frequency <= TODO
                    major_allele_freqs_by_sample = "\t".join(
                        format(-1.0 if allele == 'N' else (freq if allele == major_allele else 1.0 - freq), DECIMALS)
                        for allele, freq in site_info[9:])
                    out_sites.write(major_allele + "\t" + minor_allele + "\t" + major_allele_freqs_by_sample + "\n")


def species_worker(species_id, species_samples, args):
# def per_species(sample_names, thread_id)
    # XSNP's way: thread_id decides which banded species we would process...

    # One thing I remember about the banded: acc will run out of memory if I don't band the giant pileup files into segments ...

    outdir = f"{args.outdir}/merged/snps"
    paramstr = f"db{args.site_depth}.sc{args.sample_counts}"
    tempdir = f"{outdir}/{species_id}/temp_{paramstr}"
    if not os.path.isdir(tempdir):
        os.mkdir(tempdir)

    sample_names = species_samples[species_id]["sample_names"]

    accumulator = defaultdict(dict)
    for sample_index, sample_name  in enumerate(sample_names):
        curr_genome_coverage = species_samples[species_id]["genome_coverage"][sample_index]
        ps_args = (speices_id, sample_index, sample_name, len(sample_names), curr_genome_coverage)
        accumulate(accumulator, ps_args)

    write(accumulator, sample_names, outdir)


def midas_merge_snps(args):

    outdir = f"{args.outdir}/merged/snps"
    if args.debug and os.path.exists(outdir):
        tsprint(f"INFO:  Reusing existing output data in {outdir} according to --debug flag.")
    else:
        command(f"rm -rf {outdir}")
        command(f"mkdir -p {outdir}")

    # Read in sample list
    samples = read_samples(args.sample_list)
    # Select species
    species_samples = select_species(samples, args)

    # iggdb species_info => species_id	species_alt_id	species_name	representative_genome
    # iggdb genomes info => genome_id	genome_name	species_id
    # I want to leave create output directory not inside the select_species ...

    mp = multiprocessing.Pool(num_physical_cores)
    results = mp.map(species_worker, [(species_samples, thread_id) for thread_id in range(num_physical_cores)])

    assert all(s == "worked" for s in results), f"error"

    # About parallel computing:
    # MIDAS's way: seems to batch_samples
    # XSNP's way: each process handle banded-species
    # IGGtools's way: will we support AWS batch?
    # How about nextflow?

    # There are three tasks:
    # 1. port midas_merge_snps sub functions
    # 2. parallele compute: for each species, multiple samples and also, multiple species_samples
    # 3. possibly run out of memory issues for certain species

@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_merge_snps(args)
