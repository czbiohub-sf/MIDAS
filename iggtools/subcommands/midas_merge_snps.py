import json
import os
from collections import defaultdict
from operator import itemgetter
import multiprocessing
from math import ceil
import Bio.SeqIO

from iggtools.models.pool import Pool, select_species, search_species
from iggtools.common.utils import tsprint, num_physical_cores, command, split, InputStream, OutputStream, multiprocessing_map, multithreading_map, download_reference, select_from_tsv
from iggtools.params import outputs
from iggtools.models.uhgg import UHGG, imported_genome_file
from iggtools.params.schemas import snps_profile_schema, snps_pileup_schema, snps_info_schema, DECIMALS
from iggtools.subcommands.midas_run_snps import cat_files, scan_contigs
from iggtools.common.argparser import add_subcommand


DEFAULT_SAMPLE_COUNTS = 2
DEFAULT_GENOME_DEPTH = 5.0
DEFAULT_GENOME_COVERAGE = 0.4
DEFAULT_CHUNK_SIZE = 10000

DEFAULT_SITE_DEPTH = 1
DEFAULT_SITE_RATIO = 2.0

DEFAULT_SITE_PREV = 0.80
DEFAULT_SITE_TYPE = "common"

DEFAULT_SNP_POOLED_METHOD = "prevalence"
DEFAULT_SNP_MAF = 0.05
DEFAULT_SNP_TYPE = "mono, bi"


def create_lookup_table(packed_args):
    """ Assign or chunk up (large) contigs into chunks_of_sites for one species"""
    species_id, genome_id, chunk_size = packed_args

    global pool_of_samples

    # TODO: create the files at one place
    species_tempdir = pool_of_samples.create_species_tempdir(species_id)
    proc_lookup_file = pool_of_samples.create_species_lookup_table(species_id)

    # Stream from S3 to extract the contig_ids
    contigs_file = download_reference(imported_genome_file(genome_id, species_id, "fna.lz4"), f"{species_tempdir}")

    slice_id = 1
    with OutputStream(proc_lookup_file) as ostream, InputStream(contigs_file) as stream:
        for rec in Bio.SeqIO.parse(stream, 'fasta'):
            contig_id = rec.id
            contig_len = len(rec.seq)

            if contig_len <= chunk_size:
                # Pileup is 1-based index
                # what is the downstream used close/open situation?
                # I would prefer the same habit with python
                ostream.write("\t".join(map(str, (species_id, slice_id, contig_id, 1, contig_len))) + "\n")
                slice_id += 1
            else:
                chunk_num = ceil(contig_len/chunk_size) - 1
                for ni, ci in enumerate(range(0, contig_len, chunk_size)):
                    if ni == chunk_num:
                        ostream.write("\t".join(map(str, (species_id, slice_id, contig_id, ci+1, contig_len))) + "\n")
                    else:
                        ostream.write("\t".join(map(str, (species_id, slice_id, contig_id, ci+1, ci+chunk_size))) + "\n")
                    slice_id += 1
    return proc_lookup_file


def acgt_string(A, C, G, T):
    return ','.join(map(str, (A, C, G, T)))


def accumulate(accumulator, ps_args):
    """ Accumulate read_counts and sample_counts for a chunk of sites for one sample,
    at the same time remember <site, sample>'s A, C, G, T read counts."""

    contig_id, contig_start, contig_end, sample_index, snps_pileup_path, total_samples_count, genome_coverage = ps_args

    global global_args
    args = global_args

    # Output column indices
    c_A, c_C, c_G, c_T, c_count_samples, c_scA, c_scC, c_scG, c_scT = range(9)

    # Alternative way is to read once to memory
    awk_command = f"awk \'$1 == \"{contig_id}\" && $2 >= {contig_start} && $2 <= {contig_end}\'"
    with InputStream(snps_pileup_path, awk_command) as stream:
        for row in select_from_tsv(stream, schema=snps_pileup_schema, result_structure=dict):
            # Unpack frequently accessed columns
            ref_id, ref_pos, ref_allele = row["ref_id"], row["ref_pos"], row["ref_allele"]
            A, C, G, T, depth = row["count_a"], row["count_c"], row["count_g"], row["count_t"], row["depth"]

            # Per sample site filters:
            # if the given <site.i, sample.j> fails the within-sample site filter,
            # then sample.j should not be used for the calculation of site.i pooled statistics.
            site_ratio = depth / genome_coverage
            if depth < args.site_depth:
                continue
            if site_ratio > args.site_ratio:
                continue

            # Compute derived columns
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

            # This just remember the value from each sample.
            # Under sparse mode, site with zero read counts are not kept.
            acgt_str = acgt_string(A, C, G, T)
            assert acc[9 + sample_index] == '0,0,0,0' and acgt_str != '0,0,0,0', f"accumulate error::{site_id}:{acc}:{sample_index}"
            acc[9 + sample_index] = acgt_str


def call_alleles(alleles_all, site_depth, snp_maf):
    """ Compute the pooled allele frequencies and call SNPs """
    # Once you see all the samples, you can call SNPs
    # keep alleles passing the min allele frequency
    alleles_above_cutoff = tuple(al for al in alleles_all if al[1] / site_depth >= snp_maf)

    # classify SNPs type
    number_alleles = len(alleles_above_cutoff)
    snp_type = ["mono", "bi", "tri", "quad"][number_alleles - 1]

    # In the event of a tie -- biallelic site with 50/50 freq split -- the allele declared major is
    # the one that comes later in the "ACGT" lexicographic order.
    alleles_above_cutoff = sorted(alleles_above_cutoff, key=itemgetter(1), reverse=True)[:2]
    major_allele = alleles_above_cutoff[0][0]
    minor_allele = alleles_above_cutoff[-1][0] # for fixed sites, same as major allele

    return (major_allele, minor_allele, snp_type)


def pool_snps_and_write(accumulator, total_samples_count, species_id, slice_id, snps_info_fp, snps_freq_fp, snps_depth_fp):
    """ For each site, compute the pooled-major-alleles, site_depth, and vector of sample_depths and sample_minor_allele_freq"""

    global global_args
    args = global_args

    # Here is another places where we should created the output and temp directories
    #snps_info_fp, snps_freq_fp, snps_depth_fp = pool_of_samples.create_species_pool_snps_chunk(species_id, slice_id)

    with OutputStream(snps_info_fp) as out_info, \
            OutputStream(snps_freq_fp) as out_freq, \
                OutputStream(snps_depth_fp) as out_depth:

        for site_id, site_info in accumulator.items():
            rcA, rcC, rcG, rcT, count_samples, scA, scC, scG, scT = site_info[:9]

            # Skip site with low prevalence for core sites and vice versa for rare sites
            prevalence = count_samples / total_samples_count
            if args.site_type == "common" and prevalence < args.site_prev:
                continue
            if args.site_type == "rare" and prevalence > args.site_prev:
                continue

            # compute the pooled major allele based on the pooled-read-counts (abundance) or pooled-sample-counts (prevalence)
            if args.snp_pooled_method == "abundance":
                site_depth = sum((rcA, rcC, rcG, rcT))
                alleles_all = (('A', rcA), ('C', rcC), ('G', rcG), ('T', rcT))
            else:
                site_depth = count_samples
                alleles_all = (('A', scA), ('C', scC), ('G', scG), ('T', scT))

            major_allele, minor_allele, snp_type = call_alleles(alleles_all, site_depth, args.snp_maf)
            major_index = 'ACGT'.index(major_allele)
            minor_index = 'ACGT'.index(minor_allele)

            # Keep sites with desired snp_type
            if ('any' not in args.snp_type and snp_type not in args.snp_type):
                continue

            # Extract the read counts of pooled major alleles for samples
            sample_depths = [] # only accounts for reads matching either major or minor allele
            sample_mafs = [] # frequency of minor allele frequency
            for sample_index in range(9, len(site_info)):
                # for each <site, sample> pair
                rc_ACGT = [int(rc) for rc in site_info[sample_index].split(",")]
                if major_index == minor_index:
                    sample_depth = rc_ACGT[major_index]
                    sample_mafs.append(0.0)
                else:
                    sample_depth = rc_ACGT[major_index] + rc_ACGT[minor_index]
                    sample_mafs.append(rc_ACGT[minor_index] / sample_depth)
                sample_depths.append(sample_depth)

            # Write
            out_info.write(f"{site_id}\t{major_allele}\t{minor_allele}\t{count_samples}\t{snp_type}\t{rcA}\t{rcC}\t{rcG}\t{rcT}\t{scA}\t{scC}\t{scG}\t{scT}\n")

            write_mafs = "\t".join((str(format(maf, DECIMALS)) for maf in sample_mafs))
            out_freq.write(f"{site_id}\t" + write_mafs + "\n")

            write_depths = "\t".join((str(depth) for depth in sample_depths))
            out_depth.write(f"{site_id}\t" + write_depths + "\n")

    return True


def process_chunk(packed_args):

    global semaphore_for_species
    global slice_counts
    global list_of_species
    global species_sliced_pileup_path

    if packed_args[1] == -1:
        # return a status flag
        # the path should be computable somewhere else

        species_id = packed_args[0]
        for _ in range(slice_counts[species_id]):
            semaphore_for_species[species_id].acquire()
        print(f"++++++++++++++++++++ start {species_id}")

        # For the given species_id, what are the files I need to merge?

        samples_name = search_species(list_of_species, species_id).fetch_samples_name()
        current_species_chunks = species_sliced_pileup_path[species_id]
        print(f"merge for {species_id} => {current_species_chunks}")
        return merge_chunks_by_species(current_species_chunks, species_id, samples_name)

    try:
        # For each process, scan over all the samples
        species_id, slice_id, contig_id, contig_start, contig_end, samples_depth, samples_snps_pileup, total_samples_count = packed_args

        accumulator = dict()
        for sample_index in range(total_samples_count):
            genome_coverage = samples_depth[sample_index]
            snps_pileup_dir = samples_snps_pileup[sample_index]
            ps_args = (contig_id, contig_start, contig_end, sample_index, snps_pileup_dir, total_samples_count, genome_coverage)
            accumulate(accumulator, ps_args)

        # Compute and write pooled SNPs for each chunk of genomic sites
        snps_info_fp, snps_freq_fp, snps_depth_fp = species_sliced_pileup_path[species_id][slice_id]
        pool_snps_and_write(accumulator, total_samples_count, species_id, slice_id, snps_info_fp, snps_freq_fp, snps_depth_fp)
        #tsprint(f"process_chunk::Finished processing Species {species_id} - Process ID {slice_id} for contig {contig_id}.")
        #return {f"{species_id}_{slice_id}": (snps_info_fp, snps_freq_fp, snps_depth_fp)}
        return True

    finally:
        semaphore_for_species[species_id].release() # no deadlock
        print(f"{species_id} release once")


def merge_chunks_by_species(current_species_chunks, species_id, samples_name):

    global pool_of_samples
    species_outdir = pool_of_samples.create_species_outdir(species_id)


    # Previously generated local temp chunk files per species_id
    local_files = list(current_species_chunks.values())
    snps_info_files = [fi[0] for fi in local_files]
    snps_freq_files = [fi[1] for fi in local_files]
    snps_depth_files = [fi[2] for fi in local_files]

    snps_info_fp, snps_freq_fp, snps_depth_fp = pool_of_samples.create_species_pool_snps_results(species_id)
    # Add header for the merged-chunks
    with OutputStream(snps_info_fp) as out_info:
        out_info.write("\t".join(list(snps_info_schema.keys())) + "\n")
    cat_files(snps_info_files, snps_info_fp, chunk_num=20)

    with OutputStream(snps_freq_fp) as out_freq:
        out_freq.write("site_id\t" + "\t".join(samples_name) + "\n")
    cat_files(snps_freq_files, snps_freq_fp, chunk_num=10)

    with OutputStream(snps_depth_fp) as out_depth:
        out_depth.write("site_id\t" + "\t".join(samples_name) + "\n")
    cat_files(snps_depth_files, snps_depth_fp, chunk_num=10)

    return True


def midas_merge_snps(args):

    # tempdir indicates the input arguments
    global pool_of_samples
    global list_of_species

    pool_of_samples = Pool(args.samples_list, args.midas_outdir, "snps")
    list_of_species = select_species(pool_of_samples, "snps", args,)

    global global_args
    global_args = args

    # Create species-to-process lookup table for each species
    chunk_size = args.chunk_size
    local_toc = download_reference(outputs.genomes, pool_of_samples.get_target_layout("dbsdir"))
    db = UHGG(local_toc)

    argument_list = []
    for species in list_of_species:
        argument_list.append((species.id, db.fetch_representative_genome_id(species.id), chunk_size))
    lookups = multiprocessing_map(create_lookup_table, argument_list, num_procs=num_physical_cores)

    # dict of semaphore and acquire before
    global semaphore_for_species
    global slice_counts
    global species_sliced_pileup_path

    slice_counts = dict()
    semaphore_for_species = dict()
    species_sliced_pileup_path = defaultdict(dict)

    argument_list = []
    for si, species in enumerate(list_of_species):
        species_id = species.id
        total_samples_count = len(species.samples)
        samples_depth = species.samples_depth
        samples_snps_pileup = [sample.get_target_layout("snps_pileup", species_id) for sample in list(species.samples)]

        slice_counter = 0
        with open(lookups[si]) as stream:
            for line in stream:
                species_id, slice_id, contig_id, contig_start, contig_end = tuple(line.strip("\n").split("\t"))

                # Here is another places where we should created the output and temp directories
                snps_info_fp, snps_freq_fp, snps_depth_fp = pool_of_samples.create_species_pool_snps_chunk(species_id, slice_id)
                species_sliced_pileup_path[species_id][slice_id] = (snps_info_fp, snps_freq_fp, snps_depth_fp)

                # TODO: we need to think of another way to pass on the arguments list: several lists of hundreds of elements is too bad
                # can cause overhead passing arguments around
                my_args = (species_id, slice_id, contig_id, int(contig_start), int(contig_end), samples_depth, samples_snps_pileup, total_samples_count)
                argument_list.append(my_args)

                slice_counter += 1

            # Submit the merge jobs
            argument_list.append((species_id, -1))

        semaphore_for_species[species_id] = multiprocessing.Semaphore(slice_counter)
        slice_counts[species_id] = slice_counter
        for _ in range(slice_counter):
            semaphore_for_species[species_id].acquire()

    # Accumulate and compute pooled SNPs stastics by chunks and write tempdir
    chunks_files = multiprocessing_map(process_chunk, argument_list, num_procs=num_physical_cores)
    # Do I still missing the merge_snps_summary?
    # I should try to finish the merge fast and stretch and leave the rest for tomorrow.
    print(chunks_files)


def register_args(main_func):
    subparser = add_subcommand('midas_merge_snps', main_func, help='pooled-samples SNPs calling')

    subparser.add_argument('midas_outdir',
                           type=str,
                           help="""Path to directory to store results.  Subdirectory will be created for each species.""")
    subparser.add_argument('--samples_list',
                           dest='samples_list',
                           type=str,
                           required=True,
                           help=f"TSV file mapping sample name to midas_run_species.py output directories")
    subparser.add_argument('--species_list',
                           dest='species_list',
                           type=str,
                           metavar="CHAR",
                           help=f"Comma separated list of species ids")
    subparser.add_argument('--chunk_size',
                           dest='chunk_size',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_CHUNK_SIZE,
                           help=f"Number of genomic sites for the temporary chunk file  ({DEFAULT_CHUNK_SIZE})")


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
    # Species filters
    subparser.add_argument('--sample_counts',
                           dest='sample_counts', #min_samples
                           type=int,
                           metavar="INT",
                           default=DEFAULT_SAMPLE_COUNTS,
                           help=f"select species with >= MIN_SAMPLES ({DEFAULT_SAMPLE_COUNTS})")

    # Per sample site filters
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

    # Across samples site filters
    subparser.add_argument('--site_prev',
                           dest='site_prev',
                           default=DEFAULT_SITE_PREV,
                           type=float,
                           metavar="FLOAT",
                           help=f"Minimum fraction of sample where genomic site satifying the site filters ({DEFAULT_SITE_PREV})")
    subparser.add_argument('--site_type',
                           dest='site_type',
                           type=str,
                           default=DEFAULT_SITE_TYPE,
                           choices=['common', 'rare'],
                           help=f"Either core SNPs or rare SNPs ({DEFAULT_SITE_TYPE})")

    # SNPs calling
    subparser.add_argument('--snp_pooled_method',
                           dest='snp_pooled_method',
                           type=str,
                           default=DEFAULT_SNP_POOLED_METHOD,
                           choices=['prevalence', 'abundance'],
                           help=f"Method of call across-samples-pooled-SNPs based on either prevalence or abundance (Default: {DEFAULT_SNP_POOLED_METHOD}).")
    subparser.add_argument('--snp_maf',
                           dest='snp_maf',
                           type=float,
                           metavar="FLOAT",
                           default=DEFAULT_SNP_MAF,
                           help=f"Minimum pooled-minor-allele_frequency to call an allele present ({DEFAULT_SNP_MAF}), Values > 0.0 and < 0.5 are accepted.")
    subparser.add_argument('--snp_type',
                           type=str,
                           dest='snp_type',
                           default=DEFAULT_SNP_TYPE,
                           choices=['any', 'mono', 'bi', 'tri', 'quad'],
                           nargs='+',
                           help="""Specify one or more of the following:
                                    mono: keep sites with 1 allele > DEFAULT_SNP_MAF
                                    bi: keep sites with 2 alleles > DEFAULT_SNP_MAF
                                    tri: keep sites with 3 alleles > DEFAULT_SNP_MAF
                                    quad: keep sites with 4 alleles > DEFAULT_SNP_MAF
                                    any: keep sites regardless of observed alleles
                                    (Default: {DEFAULT_SNP_TYPE})""")

    return main_func


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_merge_snps(args)
