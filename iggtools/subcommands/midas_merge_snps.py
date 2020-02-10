import json
import os
import sys
from collections import defaultdict
from operator import itemgetter

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, num_physical_cores, command, InputStream, OutputStream, multithreading_hashmap, multithreading_map, download_reference, select_from_tsv
from iggtools.params import outputs
from iggtools.models.uhgg import UHGG


DEFAULT_SAMPLE_COUNTS = 2
DEFAULT_GENOME_DEPTH = 5.0
DEFAULT_GENOME_COVERAGE = 0.4

DEFAULT_SITE_DEPTH = 1
DEFAULT_SITE_RATIO = 2.0
DEFAULT_SITE_PREV = 0.80
DEFAULT_SITE_MAF = 0.01
DEFAULT_ALLELE_FREQ = 0.01

DEFAULT_SNP_MAF = 0
DEFAULT_SNP_TYPE = "bi"
DEFAULT_ALLELE_TYPE = "sample_counts"

DEBUG_MAX_LINES = 1000 * 1000

DECIMALS = ".3f"


## TODO: how about contigs_stats? shall we add it back to midas?
snps_summary_schema = {
    "species_id": str,
    "genome_length": int,
    "covered_bases": int,
    "total_depth": int,
    "aligned_reads": int,
    "mapped_reads": int,
    "fraction_covered": float,
    "mean_coverage": float,
}


snps_pileup_schema = {
    "ref_id": str,
    "ref_pos": int,
    "ref_allele": str,
    "depth": int,
    "count_a": int,
    "count_c": int,
    "count_g": int,
    "count_t": int,
}


snps_info_schema = {
    "site_id": str,
    "major_allele": str,
    "minor_allele": str,
    "sample_counts": int,
    "snp_type": str,
    "rc_A": int,
    "rc_C": int,
    "rc_G": int,
    "rc_T": int,
    "sc_A": int,
    "sc_C": int,
    "sc_G": int,
    "sc_T": int,
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

    # Species and sample filters
    subparser.add_argument('--genome_depth',
                           dest='genome_depth',
                           type=float,
                           metavar="FLOAT",
                           default=DEFAULT_GENOME_DEPTH,
                           help=f"Minimum average read depth per sample ({DEFAULT_GENOME_DEPTH})")
    subparser.add_argument('--genome_coverage',
                           dest='genome_coverage', #fract_cov
                           type=float,
                           metavar="FLOAT",
                           default=DEFAULT_GENOME_COVERAGE,
                           help=f"Fraction of reference sites covered by at least 1 read ({DEFAULT_GENOME_COVERAGE})")
    subparser.add_argument('--sample_counts',
                           dest='sample_counts', #min_samples
                           type=int,
                           metavar="INT",
                           default=DEFAULT_SAMPLE_COUNTS,
                           help=f"select species with >= MIN_SAMPLES ({DEFAULT_SAMPLE_COUNTS})")

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
    subparser.add_argument('--site_maf', # newly added
                           dest='site_maf',
                           default=DEFAULT_SITE_MAF,
                           type=float,
                           metavar="FLOAT",
                           help=f"Minimum minor allele frequency to call a pooled SNP site ({DEFAULT_SITE_MAF})")

    # SNPs calling
    subparser.add_argument('--allele_type',
                           dest='allele_type',
                           type=str,
                           default=DEFAULT_ALLELE_TYPE,
                           choices=['sample_counts', 'read_counts'],
                           help=f"Methods to call pooled major/minor alleles based on (Default method: {DEFAULT_ALLELE_TYPE}).")
    subparser.add_argument('--allele_freq',
                           dest='allele_freq',
                           type=float,
                           metavar="FLOAT",
                           default=DEFAULT_ALLELE_FREQ,
                           help=f"Minimum frequency for calling an allele present ({DEFAULT_ALLELE_FREQ}).")
    subparser.add_argument('--snp_maf',
                           dest='snp_maf',
                           type=float,
                           metavar="FLOAT",
                           default=DEFAULT_SNP_MAF,
                           help=f"Minimum frequency for calling an across samples pooled allele present ({DEFAULT_SNP_MAF}), Values > 0.0 and < 0.5 are accepted.")
    subparser.add_argument('--snp_type',
                           type=str,
                           dest='snp_type',
                           default=DEFAULT_SNP_TYPE,
                           choices=['any', 'mono', 'bi', 'tri', 'quad'],
                           nargs='+',
                           help="""Specify one or more of the following:
                                    mono: keep sites with 1 allele > DEFAULT_ALLELE_FREQ
                                    bi: keep sites with 2 alleles > DEFAULT_ALLELE_FREQ (default)
                                    tri: keep sites with 3 alleles > DEFAULT_ALLELE_FREQ
                                    quad: keep sites with 4 alleles > DEFAULT_ALLELE_FREQ
                                    any: keep sites regardless of observed alleles""")

    return main_func


def imported_genome_file(genome_id, species_id, component):
    return f"{outputs.cleaned_imports}/{species_id}/{genome_id}/{genome_id}.{component}"


def read_samples(sample_list):
    """
    Input: read in Table-of-Content: sample_name\tpath/to/midas/snps_dir

    Output: return dict of samples: local file path of midas_snps_summary
    """

    with InputStream(sample_list) as stream:
        samples = dict(select_from_tsv(stream, selected_columns=["sample_name", "midas_output_path"]))

    samples_midas = defaultdict()
    for sample_name in samples.keys():
        midas_output_path = samples[sample_name]
        assert os.path.exists(midas_output_path), f"MIDAS output directory {midas_output_path} for sample {sample_name} not exist."

        snps_summary = f"{midas_output_path}/snps/output/summary.txt"
        assert os.path.exists(snps_summary), f"Missing MIDAS snps profile: {snps_summary}"

        samples_midas[sample_name] = snps_summary

    return samples_midas


def read_snps_summary(snps_summary_path):
    with InputStream(snps_summary_path) as instream:
        for record in select_from_tsv(instream, selected_columns=snps_summary_schema, result_structure=dict):
            yield record


def sort_species(species_samples, sort_by="sample_counts"):
    """
    Sort species by sample_counts in descending order

    Output: list of sorted species_ids
    """

    species_ids = tuple(species_samples.keys())
    species_values = (species_samples[species_id][sort_by] for species_id in species_ids)

    sorted_species = sorted(zip(species_ids, species_values), key=itemgetter(1), reverse=True)

    return (_[0] for _ in sorted_species)


def select_species(samples, args, outdir):
    """
    Select high quality sample-species pairs based on the genome_coverage and
    genome_depth, for a list of samples; further on filter out low prevalent species
    using sample_counts cutoff.

    Input:
        - MIDAS run snps summary

    Ouput: <species, samples> pair
        - for each species, return a dictionary with two separate list: sample_name, and genome_coverage.
        - write snps_summary to file
    """

    species_samples = defaultdict(lambda: defaultdict(list))
    species_snps_summary = defaultdict(list)

    for sample_name, snps_summary_path in samples.items():
        # Read in the SNPs summary files
        for record in read_snps_summary(snps_summary_path):
            species_id = record["species_id"]
            mean_coverage = record["mean_coverage"]

            # For each <sample, species> pairs:
            # filter out low abundant species
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
            species_snps_summary[species_id].append([sample_name] + list(record.values())[1:])

    # Filter out low prevalent species
    sorted_species_ids = sort_species(species_samples, "sample_counts")
    species_samples_pass = {species_id: val for species_id, val in species_samples.items() if species_samples[species_id]["sample_counts"] >= args.sample_counts}

    # Write merged snps_summary
    pass_species_ids = species_samples_pass.keys()
    snps_summary_dir = f"{outdir}/snps_summary.tsv.lz4"

    with OutputStream(snps_summary_dir) as stream:
        stream.write("\t".join(["species_id", "sample_id"] + list(snps_summary_schema.keys())[1:]) + "\n")
        for species_id, species_list in species_snps_summary.items():
            for i in range(len(species_list)):
                stream.write("\t".join(map(str, [species_id] + species_list[i])) + "\n")

    return species_samples_pass


def read_contig_ids(contigs_file):
    contig_ids = []
    with InputStream(contigs_file) as stream:
        for line in stream:
            if line.startswith(">"):
                contig_ids.append(line.rstrip("\n")[1:])
    return contig_ids


def collect_contigs(species_samples, tempdir):
    "Generate the mapping between contigs_id - (genome_name) - species_id"

    # TODO: this should be pre-computed during the building representative genomes step.

    local_toc = download_reference(outputs.genomes)
    db = UHGG(local_toc)
    representatives = db.representatives

    def download_contigs(species_id):
        genome_id = representatives[species_id]
        return download_reference(imported_genome_file(genome_id, species_id, "fna.lz4"), f"{tempdir}/{species_id}")

    # species_id: local path to the downloaded representative genome fasta file
    species_ids = list(species_samples.keys())
    genome_files = multithreading_hashmap(download_contigs, species_ids, num_threads=20)

    # species_id: list of contig_ids
    contig_ids = multithreading_map(read_contig_ids, list(genome_files.values()))
    species_contigs = dict(zip(species_ids, contig_ids))

    return species_contigs


def acgt_string(A, C, G, T):
    return ','.join(map(str, (A, C, G, T)))


def read_snps_pileup(snps_pileup_path, contig_ids):
    with InputStream(snps_pileup_path) as instream:
        for row in select_from_tsv(instream, selected_columns=snps_pileup_schema, result_structure=dict):
            if row["ref_id"] in contig_ids:
                yield row


def accumulate(accumulator, ps_args):
    """
    Input:
        species_id: compute for ALL non-zero depth sites
        contig_ids: the list of contigs belonging to the given species
        sample_index: column index for the acuumulator matrix
        sample_name:
        snps_dir: path to MIDAS's run snps wokflow.
        total_sample_counts: total number of samples for the given species. (self.total_counts)
        genome_coverage: the precomputed mean genome coverage for <species_id, samples_list> pair.

    Output:
        compuate site_counts (rc_A, rc_C, rc_G, rc_T)

    Bascially, we want to slice by samples, instead of species_id. And for the pooled statistics,
    we need to see ALL samples before any compute.
    """

    args, contig_ids, sample_index, snps_pileup_dir, total_sample_counts, genome_coverage = ps_args

    table_iterator = read_snps_pileup(snps_pileup_dir, contig_ids)

    # Output column indices
    c_A, c_C, c_G, c_T, c_count_samples, c_scA, c_scC, c_scG, c_scT = range(9)

    for line, row in enumerate(table_iterator):
        if line == DEBUG_MAX_LINES:
            break

        # Unpack frequently accessed columns
        contig_id = row["ref_id"]
        ref_pos = row["ref_pos"]
        ref_allele = row["ref_allele"]
        depth = row["depth"]
        A, C, G, T = row["count_a"], row["count_c"], row["count_g"], row["count_t"]

        # Per sample site filters:
        # if the given <site.i, sample.j> fails the within-sample site filter,
        # then sample.j should not be used for the calculation of site.i pooled statistics.
        site_ratio = depth / genome_coverage
        if depth < args.site_depth:
            continue
        if site_ratio > args.site_ratio:
            continue

        # Computae derived columns
        site_id = f"{contig_id}|{ref_pos}|{ref_allele}"

        # sample counts for A, C, G, T
        sc_ACGT = [0, 0, 0, 0]
        for i, nt_count in enumerate((A, C, G, T)):
            if nt_count > 0: # presence or absence
                sc_ACGT[i] = 1

        # Aggragate
        acc = accumulator.get(site_id)
        if acc:
            acc[c_A] += A
            acc[c_C] += C
            acc[c_G] += G
            acc[c_T] += T
            acc[c_count_samples] += 1
            acc[c_scA] += sc_ACGT[0]
            acc[c_scC] += sc_ACGT[1]
            acc[c_scG] += sc_ACGT[2]
            acc[c_scT] += sc_ACGT[3]
        else:
            # initialize each sample_index column with 0,0,0,0, particularly
            # for <site, sample> pair either absent or fail the site filters
            acc = [A, C, G, T, 1, sc_ACGT[0], sc_ACGT[1], sc_ACGT[2], sc_ACGT[3]] + ([acgt_string(0, 0, 0, 0)] * total_sample_counts)
            accumulator[site_id] = acc

        # This isn't being accumulated across samples;
        # we are just remembering the value from each sample.
        acgt_str = acgt_string(A, C, G, T)
        assert acc[9 + sample_index] == '0,0,0,0' and acgt_str != '0,0,0,0'
        acc[9 + sample_index] = acgt_str


def call_alleles(alleles_all, site_depth, site_maf):
    """
    After seeing data across ALL the samples, we can compute the pooled allele frequencies and call alleles.
    """
    # keep alleles passing the min allele frequency
    alleles_above_cutoff = tuple(al for al in alleles_all if al[1] / site_depth >= site_maf)

    # classify SNPs type
    number_alleles = len(alleles_above_cutoff)
    snp_type = ["mono", "bi", "tri", "quad"][number_alleles - 1]

    if number_alleles > 0:
        # In the event of a tie -- biallelic site with 50/50 freq split -- the allele declared major is
        # the one that comes later in the "ACGT" lexicographic order.
        alleles_above_cutoff = sorted(alleles_above_cutoff, reverse=True)[:2]
        major_allele = alleles_above_cutoff[0][0]
        minor_allele = alleles_above_cutoff[-1][0] # for fixed sites, same as major allele
    else:
        sys.exit("Error: no alleles for pooled site: This should not happen")

    return (major_allele, minor_allele, snp_type)


def pool_and_write(accumulator, sample_names, outdir, args):
    """
    Ask Boris: write intermediate accumulator to file, and then re-read in here?
        - I/O bandwidth
        + no need to re-accumulate with potentially downstream bugs

    Output: after seeing ALL samples' pileup data, we compute the pooled statistics for each site
        - start with pooled_major_alleles and pooled_minor_alleles
        - then easily site_counts and site_depth
        - then calcualte the per sample list for sample_depths
        - then sample_mafs
    """

    snps_info_fp = f"{outdir}/snps_info.tsv.lz4"
    snps_freq_fp = f"{outdir}/snps_freqs.tsv.lz4"
    snps_depth_fp = f"{outdir}/snps_depth.tsv.lz4"

    with OutputStream(snps_freq_fp) as stream_freq, \
            OutputStream(snps_depth_fp) as stream_depth, \
                OutputStream(snps_info_fp) as stream_info:

        stream_info.write("\t".join(list(snps_info_schema.keys())) + "\n")
        stream_freq.write("site_id\t" + "\t".join(sample_names) + "\n")
        stream_depth.write("site_id\t" + "\t".join(sample_names) + "\n")

        for site_id, site_info in accumulator.items():

            rcA, rcC, rcG, rcT, count_samples, scA, scC, scG, scT = site_info[:9]

            # skip site with low prevalence
            prevalence = count_samples / len(sample_names)
            if prevalence < args.site_prev:
                continue

            # compute the pooled major allele based on the pooled-read-counts (abundance) or pooled-sample-counts (prevalence)
            if args.allele_type == "read_counts":
                site_counts = (rcA, rcC, rcG, rcT)
                alleles_all = (('A', rcA), ('C', rcC), ('G', rcG), ('T', rcT))
            else:
                site_counts = (scA, scC, scG, scT)
                alleles_all = (('A', scA), ('C', scC), ('G', scG), ('T', scT))

            site_depth = sum(site_counts)
            major_allele, minor_allele, snp_type = call_alleles(alleles_all, site_depth, args.site_maf)
            major_index = 'ACGT'.index(major_allele)
            minor_index = 'ACGT'.index(minor_allele)

            # TODO: need to handle the 'any' as snp_types input arguments
            snp_types = ['quad', 'tri', 'bi', 'mono']
            if ('any' not in snp_types and snp_type not in snp_types):
                continue

            sample_depths = [] # only accounts for reads matching either major or minor allele
            sample_mafs = [] # frequency of minor allele frequency
            print(len(site_info))
            for sample_index in range(9, len(site_info)):
                # for <site, sample>
                print(site_info[sample_index])
                rc_ACGT = [int(rc) for rc in site_info[sample_index].split(",")]
                sample_depths.append(rc_ACGT[major_index] + rc_ACGT[minor_index])
                sample_mafs.append(rc_ACGT[minor_index])

            print(sample_depths)
            print(sample_mafs)
            # write
            stream_info.write(f"{site_id}\t{major_allele}\t{minor_allele}\t{count_samples}\t{snp_type}\t{rcA}\t{rcC}\t{rcG}\t{rcT}\t{scA}\t{scC}\t{scG}\t{scT}\n")

            write_mafs = "\t".join((str(format(maf, DECIMALS)) for maf in sample_mafs))
            stream_freq.write(f"{site_id}\t" + write_mafs + "\n")

            write_depths = "\t".join((str(format(depth, DECIMALS)) for depth in sample_depths))
            stream_depth.write(f"{site_id}\t" + write_depths + "\n")

    return "done"


def per_species_worker(arguments):

    species_id, species_sample, contig_ids, samples_midas, args, outdir = arguments

    sample_names = species_sample["sample_names"]
    genome_coverage = species_sample["genome_coverage"]

    outdir = f"{outdir}/{species_id}"
    command(f"rm -rf {outdir}")
    command(f"mkdir -p {outdir}")

    # Scan over ALL the samples passing the genome filters
    accumulator = dict()
    for sample_index, sample_name  in enumerate(sample_names):
        sample_genome_cov = genome_coverage[sample_index]

        snps_summary_dir = samples_midas[sample_name]
        snps_dir = os.path.dirname(snps_summary_dir)
        snps_pileup_dir = os.path.join(snps_dir, f"{species_id}.snps.lz4")
        assert os.path.exists(snps_pileup_dir)

        total_sample_counts = len(sample_names)
        ps_args = (args, contig_ids, sample_index, snps_pileup_dir, total_sample_counts, sample_genome_cov)

        accumulate(accumulator, ps_args)

    flag = pool_and_write(accumulator, sample_names, outdir, args)
    assert flag == "done"


def midas_merge_snps(args):

    outdir = f"{args.outdir}/merged/snps"
    if args.debug and os.path.exists(outdir):
        tsprint(f"INFO:  Reusing existing output data in {outdir} according to --debug flag.")
    else:
        command(f"rm -rf {outdir}")
        command(f"mkdir -p {outdir}")


    paramstr = f"sd{args.site_depth}.sr{args.site_ratio}.sp{args.site_prev}.sm{args.snp_maf}"
    tempdir = f"{outdir}/temp_{paramstr}"
    if args.debug and os.path.exists(tempdir):
        tsprint(f"INFO:  Reusing existing temp intermediate data in {tempdir} according to --debug flag.")
    else:
        command(f"rm -rf {tempdir}")
        command(f"mkdir -p {tempdir}")


    # Read in table of content: list of samples
    samples_midas = read_samples(args.sample_list)

    # Select species
    species_samples = select_species(samples_midas, args, outdir)

    # Collect contig_ids for each species
    species_contigs = collect_contigs(species_samples, outdir)

    # Accumulate read_counts and sample_counts across ALL the sample passing the genome filters;
    # and at the same time remember <site, sample>'s A, C, G, T read counts.
    argument_list = []
    for species_id in list(species_samples.keys())[:1]:
        argument_list.append((species_id, species_samples[species_id], species_contigs[species_id], samples_midas, args, outdir))

    multithreading_map(per_species_worker, argument_list, num_threads=num_physical_cores)


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_merge_snps(args)
