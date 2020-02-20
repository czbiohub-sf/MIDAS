import json
import os
import sys
from collections import defaultdict
from operator import itemgetter
import Bio.SeqIO
from math import ceil

from iggtools.models.pool import Pool, select_species

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, num_physical_cores, command, InputStream, OutputStream, retry, multiprocessing_map, multithreading_map, download_reference, select_from_tsv, TimedSection
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


# The Bio.SeqIO.parse() code is CPU-bound and thus it's best to run this function in a separate process for every genome.
def create_lookup_table(packed_args):
    """ Assign or chunk up (large) contigs into chunks for one species"""
    species_id, genome_id, chunk_size, tempdir = packed_args
    if not os.path.exists(f"{tempdir}/{species_id}"):
        command(f"mkdir -p {tempdir}/{species_id}")

    # Stream from S3 to extract the contig_ids
    contigs_file = download_reference(imported_genome_file(genome_id, species_id, "fna.lz4"), f"{tempdir}/{species_id}")
    proc_lookup_file = f"{tempdir}/{species_id}/procid_lookup.tsv"

    proc_id = 1
    with OutputStream(proc_lookup_file) as ostream, InputStream(contigs_file) as stream:
        for rec in Bio.SeqIO.parse(stream, 'fasta'):
            contig_id = rec.id
            contig_len = len(rec.seq)

            if contig_len <= chunk_size:
                # Pileup is 1-based index
                ostream.write("\t".join(map(str, (proc_id, contig_id, 1, contig_len))) + "\n")
                proc_id += 1
            else:
                chunk_num = ceil(contig_len/chunk_size) - 1
                for ni, ci in enumerate(range(0, contig_len, chunk_size)):
                    if ni == chunk_num:
                        ostream.write("\t".join(map(str, (proc_id, contig_id, ci+1, contig_len))) + "\n")
                    else:
                        ostream.write("\t".join(map(str, (proc_id, contig_id, ci+1, ci + chunk_size))) + "\n")
                    proc_id += 1
    return proc_lookup_file


def acgt_string(A, C, G, T):
    return ','.join(map(str, (A, C, G, T)))


def accumulate(accumulator, ps_args):
    """ Accumulate read_counts and sample_counts for a chunk of sites for one sample"""

    contig_id, contig_start, contig_end, sample_index, snps_pileup_path, total_samples_count, genome_coverage = ps_args

    global global_args
    args = global_args

    # Output column indices
    c_A, c_C, c_G, c_T, c_count_samples, c_scA, c_scC, c_scG, c_scT = range(9)

    # Alternative wayis to read once to memory
    awk_command = f"awk \'$1 == \"{contig_id}\" && $2 >= {contig_start} && $2 < {contig_end}\'"
    with InputStream(snps_pileup_path, awk_command) as stream:
        for row in select_from_tsv(stream, schema=snps_pileup_schema, result_structure=dict):
            # Unpack frequently accessed columns
            ref_id, ref_pos, ref_allele = row["ref_id"], row["ref_pos"], row["ref_allele"]
            A, C, G, T = row["count_a"], row["count_c"], row["count_g"], row["count_t"]
            depth = row["depth"]

            # Only process current chunk of sites
            if ref_id == contig_id and ref_pos >= contig_start and ref_pos < contig_end:

                # Per sample site filters:
                # if the given <site.i, sample.j> fails the within-sample site filter,
                # then sample.j should not be used for the calculation of site.i pooled statistics.
                site_ratio = depth / genome_coverage
                if depth < args.site_depth:
                    continue
                if site_ratio > args.site_ratio:
                    continue

                # Computae derived columns
                site_id = f"{ref_id}|{ref_pos}|{ref_allele}"

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
                stream.ignore_errors()

            # This isn't being accumulated across samples;
            # we are just remembering the value from each sample.
            acgt_str = acgt_string(A, C, G, T)
            assert acc[9 + sample_index] == '0,0,0,0' and acgt_str != '0,0,0,0'
            acc[9 + sample_index] = acgt_str


def process_chunk(packed_args):

    species_id, proc_id, contig_id, contig_start, contig_end, samples_name, samples_depth, samples_data_dir, total_samples_count, species_outdir = packed_args

    # For each process, scan over all the samples
    accumulator = dict()
    for sample_index, sample_name in enumerate(samples_name):
        genome_coverage = samples_depth[sample_index]
        snps_pileup_dir = f"{samples_data_dir[sample_index]}/output/{species_id}.snps.lz4"
        assert os.path.exists(snps_pileup_dir), f"Missing pileup results {snps_pileup_dir} for {sample_name}"

        ps_args = (contig_id, contig_start, contig_end, sample_index, snps_pileup_dir, total_samples_count, genome_coverage)
        accumulate(accumulator, ps_args)
    print(accumulator)

    # Write pooled statistics for each chunk
    #flag = pool_and_write(accumulator, sample_names, species_outdir)
    tsprint(f"Process ID {proc_id} for species {species_id} for contig {contig_id} from {contig_start} to {contig_end} is done")
    #assert flag == "done"


def per_species_work(species, procid_lookup_table, outdir, tempdir):

    species_id = species.id

    species_outdir = f"{outdir}/{species_id}"
    command(f"rm -rf {species_outdir}")
    command(f"mkdir -p {species_outdir}")

    pool_of_samples = list(species.samples)
    samples_name = [sample.sample_name for sample in pool_of_samples]
    samples_depth = species.samples_depth
    samples_data_dir = [sample.data_dir for sample in pool_of_samples]
    total_samples_count = len(pool_of_samples)

    argument_list = []
    with open(procid_lookup_table) as stream:
        for line in stream:
            proc_id, contig_id, contig_start, contig_end = tuple(line.strip("\n").split("\t"))
            argument_list.append((species_id, proc_id, contig_id, int(contig_start), int(contig_end), samples_name, samples_depth, samples_data_dir, total_samples_count, species_outdir))

    multiprocessing_map(process_chunk, argument_list, num_procs=num_physical_cores) #multithreading_map


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


def midas_merge_snps(args):

    outdir, tempdir = check_outdir(args)
    pool_of_samples = Pool(args.sample_list, "snps")
    list_of_species = select_species(pool_of_samples, args, "snps")

    with TimedSection("Create lookup table for each species"):
        # Generate the mapping between contigs_id - [genome_name] - species_id
        chunk_size = 10000
        local_toc = download_reference(outputs.genomes)
        db = UHGG(local_toc)
        representative = db.representatives

        # Write the process_id to chunks of contigs for each species
        lookups = multiprocessing_map(create_lookup_table, ((species.id, representative[species.id], chunk_size, tempdir) for species in list_of_species), num_procs=num_physical_cores)


    # Under "sparse" mode, different samples' pileup files

    global global_args
    global_args = args

    for si, species in enumerate(list_of_species[:1]):
        tsprint(f"midas_merge_snps now processing {species.id}")
        per_species_work(species, lookups[si], outdir, tempdir)
        # We can pass on the lookups[si]
        #for sample_index, sample in enumerate(species.samples):
        #    print(sample_index, sample.sample_name)
        exit(0)

    # Accumulate read_counts and samples_count across ALL the sample passing the genome filters;
    # and at the same time remember <site, sample>'s A, C, G, T read counts.

    #multithreading_map(per_species_work, argument_list, num_threads=num_physical_cores)


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_merge_snps(args)
