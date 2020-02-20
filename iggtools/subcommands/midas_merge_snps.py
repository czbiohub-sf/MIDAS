import json
import os
import sys
from collections import defaultdict
from operator import itemgetter
import Bio.SeqIO

from iggtools.models.pool import Pool, select_species

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, num_physical_cores, command, InputStream, OutputStream, multithreading_hashmap, multithreading_map, download_reference, select_from_tsv
from iggtools.params import outputs
from iggtools.models.uhgg import UHGG
from iggtools.params.schemas import snps_profile_schema, snps_pileup_schema, snps_info_schema, DECIMALS


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
DEFAULT_ALLELE_TYPE = "samples_count"

DEBUG_MAX_LINES = 1000 * 1000


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
    subparser.add_argument('--species_id',
                           dest='species_id',
                           type=str,
                           metavar="CHAR",
                           help=f"Comma separated list of species ids")

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
                           choices=['samples_count', 'read_counts'],
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


def check_outdir(args):
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

    return (outdir, tempdir)


def imported_genome_file(genome_id, species_id, component):
    return f"{outputs.cleaned_imports}/{species_id}/{genome_id}/{genome_id}.{component}"


def acgt_string(A, C, G, T):
    return ','.join(map(str, (A, C, G, T)))


def accumulate(accumulator, ps_args):
    """
    Input:
        species_id: compute for ALL non-zero depth sites
        contig_ids: the list of contigs belonging to the given species
        sample_index: column index for the acuumulator matrix
        sample_name:
        snps_dir: path to MIDAS's run snps wokflow.
        total_samples_count: total number of samples for the given species. (self.total_counts)
        genome_coverage: the precomputed mean genome coverage for <species_id, samples_list> pair.

    Output:
        compuate site_counts (rc_A, rc_C, rc_G, rc_T)

    Bascially, we want to slice by samples, instead of species_id. And for the pooled statistics,
    we need to see ALL samples before any compute.
    """

    global global_args
    args = global_args

    contig_ids, sample_index, snps_pileup_dir, total_samples_count, genome_coverage = ps_args

    # Output column indices
    c_A, c_C, c_G, c_T, c_count_samples, c_scA, c_scC, c_scG, c_scT = range(9)

    # Alternative way (which might be better) is to read once for all
    with InputStream(snps_pileup_path) as stream:
        for row in select_from_tsv(stream, selected_columns=snps_pileup_schema, result_structure=dict):
            if row["ref_id"] not in contig_ids:
                continue

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
                acc = [A, C, G, T, 1, sc_ACGT[0], sc_ACGT[1], sc_ACGT[2], sc_ACGT[3]] + ([acgt_string(0, 0, 0, 0)] * total_samples_count)
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

    # In the event of a tie -- biallelic site with 50/50 freq split -- the allele declared major is
    # the one that comes later in the "ACGT" lexicographic order.
    alleles_above_cutoff = sorted(alleles_above_cutoff, key=itemgetter(1), reverse=True)[:2]
    major_allele = alleles_above_cutoff[0][0]
    minor_allele = alleles_above_cutoff[-1][0] # for fixed sites, same as major allele

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

            # Keep sites with desired snp_type
            if ('any' not in args.snp_type and snp_type not in args.snp_type):
                continue

            # Extract the read counts of pooled major alleles for samples
            sample_depths = [] # only accounts for reads matching either major or minor allele
            sample_mafs = [] # frequency of minor allele frequency
            for sample_index in range(9, len(site_info)):
                # for <site, sample> pair
                rc_ACGT = [int(rc) for rc in site_info[sample_index].split(",")]
                sample_depth = rc_ACGT[major_index] + rc_ACGT[minor_index]
                sample_depths.append(sample_depth)
                sample_mafs.append(rc_ACGT[minor_index] / sample_depth)

            # Write
            stream_info.write(f"{site_id}\t{major_allele}\t{minor_allele}\t{count_samples}\t{snp_type}\t{rcA}\t{rcC}\t{rcG}\t{rcT}\t{scA}\t{scC}\t{scG}\t{scT}\n")

            write_mafs = "\t".join((str(format(maf, DECIMALS)) for maf in sample_mafs))
            stream_freq.write(f"{site_id}\t" + write_mafs + "\n")

            write_depths = "\t".join((str(depth) for depth in sample_depths))
            stream_depth.write(f"{site_id}\t" + write_depths + "\n")

    return "done"


def per_species_worker(species, outdir, tempdir):

    global global_args
    global representative
    args = global_args

    # Generate the mapping between contigs_id - [genome_name] - species_id
    species_id = species.id
    genome_id = representative[species_id]

    # Stream from S3 to extract the contig_ids
    contigs_file = imported_genome_file(genome_id, species_id, "fna.lz4")
    with InputStream(contigs_file, "grep '^>'") as stream:
        contig_ids = [line.rstrip("\n")[1:] for line in stream]

    pool_of_samples = list(species.samples)
    genome_coverage = species.samples_depth

    outdir = f"{outdir}/{species_id}"
    command(f"rm -rf {outdir}")
    command(f"mkdir -p {outdir}")

    # Scan over ALL the samples passing the genome filters
    accumulator = dict()
    for sample_index, sample  in enumerate(pool_of_samples):
        sample.sample_name
        sample.midas_outdir
        pileup_dir = f"{sample.data_dir}/{species_id}.snps.lz4"
        assert os.path.exists(snps_pileup_dir), f"Missing pileup results {pileup_dir} for {sample.sample_name}"

        sample_genome_cov = species.samples_depth[sample_index]

        total_samples_count = len(sample_names)
        ps_args = (args, contig_ids, sample_index, snps_pileup_dir, total_samples_count, sample_genome_cov)
        # Don't pass on too many variables
        accumulate(accumulator, ps_args)

    flag = pool_and_write(accumulator, sample_names, outdir, args)
    assert flag == "done"


def midas_merge_snps(args):

    outdir, tempdir = check_outdir(args)
    pool_of_samples = Pool(args.sample_list, "snps")
    list_of_species = select_species(pool_of_samples, args, "snps")

    # Under "sparse" mode, different samples' pileup files
    # Collect contig_ids for each species
    #species_contigs = collect_contigs(species_samples, outdir)

    global global_args
    global_args = args
    global representative

    local_toc = download_reference(outputs.genomes)
    db = UHGG(local_toc)
    representative = db.representatives

    for species in list_of_species:
        tsprint(f"midas_merge_snps now processing {species.id}")
        for sample_index, sample in enumerate(species.samples):
            print(sample_index, sample.sample_name)
            #per_species_worker(species, outdir, tempdir)
    exit(0)


    # Accumulate read_counts and samples_count across ALL the sample passing the genome filters;
    # and at the same time remember <site, sample>'s A, C, G, T read counts.
    argument_list = []
    for species_id in list(species_samples.keys()):
        argument_list.append((species_id, species_samples[species_id], species_contigs[species_id], samples_midas, args, outdir))

    multithreading_map(per_species_worker, argument_list, num_threads=num_physical_cores)


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_merge_snps(args)
