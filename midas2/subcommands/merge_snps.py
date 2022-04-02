#!/usr/bin/env python3
import os
import json
import multiprocessing
from collections import defaultdict

from midas2.models.samplepool import SamplePool
from midas2.common.utils import tsprint, num_physical_cores, command, InputStream, OutputStream, multiprocessing_map, select_from_tsv, cat_files, multithreading_map
from midas2.common.utilities import annotate_site, acgt_string, scan_gene_feature, scan_fasta, compute_gene_boundary
from midas2.common.snvs import call_alleles
from midas2.models.midasdb import MIDAS_DB
from midas2.params.schemas import snps_pileup_schema, snps_pileup_basic_schema, snps_info_schema, format_data
from midas2.common.argparser import add_subcommand
from midas2.params.inputs import MIDASDB_NAMES
from midas2.models.species import load_chunks_cache


DEFAULT_SAMPLE_COUNTS = 2
DEFAULT_GENOME_DEPTH = 5.0
DEFAULT_GENOME_COVERAGE = 0.4
DEFAULT_CHUNK_SIZE = 1000000
DEFAULT_NUM_CORES = 16

DEFAULT_SITE_DEPTH = 2
DEFAULT_SITE_RATIO = 5.0

DEFAULT_SITE_PREV = 0.80
DEFAULT_SITE_TYPE = "common"

DEFAULT_SNP_POOLED_METHOD = "prevalence"
DEFAULT_SNP_MAF = 0.05
DEFAULT_SNP_TYPE = "bi, tri, quad"


def register_args(main_func):
    subparser = add_subcommand('merge_snps', main_func, help='pooled-samples SNPs calling')

    subparser.add_argument('midas_outdir',
                           type=str,
                           help="""Path to directory to store results.  Subdirectory will be created for each species.""")
    subparser.add_argument('--samples_list',
                           dest='samples_list',
                           type=str,
                           required=True,
                           help=f"TSV file mapping sample name to run_species.py output directories")
    subparser.add_argument('--chunk_size',
                           dest='chunk_size',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_CHUNK_SIZE,
                           help=f"Number of genomic sites for the temporary chunk file  ({DEFAULT_CHUNK_SIZE})")

    subparser.add_argument('--midasdb_name',
                           dest='midasdb_name',
                           type=str,
                           default="uhgg",
                           choices=MIDASDB_NAMES,
                           help=f"MIDAS Database name.")
    subparser.add_argument('--midasdb_dir',
                           dest='midasdb_dir',
                           type=str,
                           default="midasdb",
                           help=f"Local MIDAS Database path mirroing S3.")
    subparser.add_argument('--num_cores',
                           dest='num_cores',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_NUM_CORES,
                           help=f"Number of physical cores to use ({DEFAULT_NUM_CORES})")

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
                                    (Default: {%s})""" % DEFAULT_SNP_TYPE)

    subparser.add_argument('--advanced',
                           action='store_true',
                           default=False,
                           help=f"Report majore/minor allele for each genomic sites.")
    subparser.add_argument('--robust_chunk',
                           action='store_true',
                           default=False,
                           help=f"Adjust chunk_size based on species's prevalence.")

    return main_func


def calculate_chunk_size(samples_count, chunk_size):
    if samples_count <= 20:
        return 0
    if samples_count <= 100:
        return 2000000
    if samples_count > 800:
        return 100000
    if samples_count > 500:
        return 200000
    if samples_count > 200:
        return 500000
    return chunk_size


def in_place(species_counts):
    return species_counts < 50


def design_chunks_per_species(args):
    global global_args
    sp, midas_db = args
    samples_count = sp.samples_count

    if global_args.robust_chunk:
        chunk_size = calculate_chunk_size(samples_count, global_args.chunk_size)
    else:
        chunk_size = global_args.chunk_size

    if chunk_size == 0:
        return None # whole species, no need to design chunks

    sp.chunk_size = chunk_size
    return sp.compute_snps_chunks(midas_db, chunk_size, "merge")


def design_chunks(species_ids_of_interest, midas_db):
    global pool_of_samples
    global dict_of_species
    global global_args

    global semaphore_for_species
    semaphore_for_species = dict()

    global dict_of_site_chunks

    # Design chunks structure per species
    num_cores = min(midas_db.num_cores, 16)
    all_site_chunks = multithreading_map(design_chunks_per_species, [(sp, midas_db) for sp in dict_of_species.values()], num_cores) #<---

    if in_place(len(species_ids_of_interest)):
        dict_of_site_chunks = defaultdict(dict)
        for spidx, species_id in enumerate(species_ids_of_interest):
            if all_site_chunks[spidx] is not None:
                dict_of_site_chunks[species_id] = all_site_chunks[spidx]

    arguments_list = []
    for sp in dict_of_species.values():
        species_id = sp.id
        num_of_chunks = sp.num_of_snps_chunks

        if num_of_chunks is not None:
            for chunk_id in range(0, num_of_chunks):
                arguments_list.append((species_id, chunk_id))

            semaphore_for_species[species_id] = multiprocessing.Semaphore(num_of_chunks)
            for _ in range(num_of_chunks):
                semaphore_for_species[species_id].acquire()
        else:
            arguments_list.append((species_id, -2)) # species_worker

    tsprint("================= Total number of compute chunks: " + str(len(arguments_list)))

    for species_id in dict_of_species.keys():
        if species_id in semaphore_for_species:
            arguments_list.append((species_id, -1))
    return arguments_list


def process(packed_args):

    species_id, chunk_id = packed_args

    if chunk_id == -1:

        global semaphore_for_species
        global dict_of_species
        sp = dict_of_species[species_id]

        tsprint(f"  MIDAS2::process::{species_id}-{chunk_id}::wait collect_chunks")
        for _ in range(sp.num_of_snps_chunks):
            semaphore_for_species[species_id].acquire()
        tsprint(f"  MIDAS2::process::{species_id}-{chunk_id}::start collect_chunks")
        collect_chunks(species_id)
        tsprint(f"  MIDAS2::process::{species_id}-{chunk_id}::finish collect_chunks")
        return "worked"

    tsprint(f"  MIDAS2::process::{species_id}-{chunk_id}::start snps_worker")
    snps_worker(species_id, chunk_id)
    tsprint(f"  MIDAS2::process::{species_id}-{chunk_id}::finish snps_worker")
    return "worked"


def snps_worker(species_id, chunk_id):
    """ For genome sites from one chunk, scan across all the sample, compute pooled SNPs and write to file """

    global semaphore_for_species
    global dict_of_species
    global dict_of_site_chunks
    global global_args

    try:
        sp = dict_of_species[species_id]

        if chunk_id == -2:
            species_worker(species_id)
        else:
            if in_place(len(dict_of_species)):
                chunks_of_sites = dict_of_site_chunks[species_id]
            else:
                chunks_of_sites = load_chunks_cache(sp.chunks_of_sites_fp)
            chunk_worker(chunks_of_sites[chunk_id][0])
    finally:
        if species_id in semaphore_for_species:
            semaphore_for_species[species_id].release() # no deadlock


def species_worker(species_id):
    global dict_of_species

    sp = dict_of_species[species_id]

    total_samples_count = sp.samples_count
    list_of_samples_depth = sp.list_of_samples_depth
    list_of_samples = sp.list_of_samples

    tsprint(f"    MIDAS2::species_worker::{species_id}--2::start accumulate_samples")
    accumulator = dict()
    for sample_index, sample in enumerate(list_of_samples):
        snps_pileup_path = sample.get_target_layout("snps_pileup", species_id)
        proc_args = ("species", sample_index, snps_pileup_path, total_samples_count, list_of_samples_depth[sample_index])
        accumulate(accumulator, proc_args)
    tsprint(f"    MIDAS2::species_worker::{species_id}--2::finish accumulate_samples")

    tsprint(f"    MIDAS2::species_worker::{species_id}--2::start call_and_write_population_snps")
    pooled_snps_dict = call_population_snps(accumulator, species_id)
    write_population_snps(pooled_snps_dict, species_id, -2)
    tsprint(f"    MIDAS2::species_worker::{species_id}--2::finish call_and_write_population_snps")


def chunk_worker(packed_args):
    """ Accumulate sample by sample and filter population SNPs """
    global pool_of_samples
    global dict_of_species

    species_id, chunk_id, contig_id = packed_args[:3]

    sp = dict_of_species[species_id]
    total_samples_count = sp.samples_count
    list_of_samples_depth = sp.list_of_samples_depth

    tsprint(f"    MIDAS2::chunk_worker::{species_id}-{chunk_id}::start accumulate_samples")
    accumulator = dict()
    for sample_index, sample in enumerate(sp.list_of_samples):
        snps_pileup_path = sample.get_target_layout("snps_pileup", species_id)

        if contig_id == -1:
            loc_fp = sp.chunks_contigs[chunk_id]
            proc_args = ("file", sample_index, snps_pileup_path, total_samples_count, list_of_samples_depth[sample_index], loc_fp)
        else:
            # Pileup is 1-based index, close left close right
            contig_start, contig_end = packed_args[3:5]
            proc_args = ("range", sample_index, snps_pileup_path, total_samples_count, list_of_samples_depth[sample_index], contig_id, contig_start+1, contig_end)
        accumulate(accumulator, proc_args)
    tsprint(f"    MIDAS2::chunk_worker::{species_id}-{chunk_id}::finish accumulate_samples")

    # Compute across-samples SNPs and write to chunk file
    tsprint(f"    MIDAS2::chunk_worker::{species_id}-{chunk_id}::start call_and_write_population_snps")
    pooled_snps_dict = call_population_snps(accumulator, species_id)
    write_population_snps(pooled_snps_dict, species_id, chunk_id)
    tsprint(f"    MIDAS2::chunk_worker::{species_id}-{chunk_id}::finish call_and_write_population_snps")


def accumulate(accumulator, proc_args):
    """ Accumulate read_counts and sample_counts for a chunk of sites for one sample,
    at the same time remember <site, sample>'s A, C, G, T read counts."""

    global global_args

    flag, sample_index, snps_pileup_path, total_samples_count, genome_coverage = proc_args[:5]

    if flag == "file":
        loc_fp = proc_args[5]
        filter_cmd = f"grep -Fwf {loc_fp}"
    if flag == "range":
        contig_id, contig_start, contig_end = proc_args[5:]
        filter_cmd = f"awk \'$1 == \"{contig_id}\" && $2 >= {contig_start} && $2 <= {contig_end}\'"
    if flag == "species":
        filter_cmd = f"tail -n +2"

    # Output column indices
    c_A, c_C, c_G, c_T, c_count_samples, c_scA, c_scC, c_scG, c_scT = range(9)

    curr_schema = snps_pileup_schema if global_args.advanced else snps_pileup_basic_schema

    with InputStream(snps_pileup_path, filter_cmd) as stream:
        for row in select_from_tsv(stream, schema=curr_schema, selected_columns=snps_pileup_basic_schema, result_structure=dict):
            # Unpack frequently accessed columns
            ref_id, ref_pos, ref_allele = row["ref_id"], row["ref_pos"], row["ref_allele"]
            A, C, G, T, depth = row["count_a"], row["count_c"], row["count_g"], row["count_t"], row["depth"]

            # Per Sample Site Filters: if the given <site.i, sample.j> pair fails the within-sample site filter,
            # then sample.j should not be used for the calculation of site.i pooled statistics.
            site_ratio = depth / genome_coverage

            if depth < global_args.site_depth:
                continue
            if site_ratio > global_args.site_ratio:
                continue

            # Compute derived columns
            site_id = f"{ref_id}|{ref_pos}|{ref_allele}"

            # Sample counts for A, C, G, T
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
                # Initialize each sample_index column with 0,0,0,0, particularly
                # for <site, sample> pair either absent or fail the site filters
                acc = [A, C, G, T, 1, sc_ACGT[0], sc_ACGT[1], sc_ACGT[2], sc_ACGT[3]] + ([acgt_string(0, 0, 0, 0)] * total_samples_count)
                accumulator[site_id] = acc

            # This just remember the value from each sample.
            # Under sparse mode, site with zero read counts are not kept.
            acgt_str = acgt_string(A, C, G, T)
            assert acc[9 + sample_index] == '0,0,0,0' and acgt_str != '0,0,0,0', f"accumulate error::{site_id}:{acc}:{sample_index}:{acgt_str}"
            acc[9 + sample_index] = acgt_str

        stream.ignore_errors()


def call_population_snps(accumulator, species_id):
    """ For each site, compute the pooled-major-alleles, site_depth, and vector of sample_depths and sample_minor_allele_freq"""

    global global_args
    global dict_of_species

    sp = dict_of_species[species_id]
    total_samples_count = sp.samples_count

    genes_feature = scan_gene_feature(sp.gene_feature_fp)
    genes_sequence = scan_fasta(sp.gene_seq_fp)
    genes_boundary = compute_gene_boundary(genes_feature)

    pooled_snps_dict = {
        "info": {},
        "freq": {},
        "depth": {},
    }

    for site_id, site_info in accumulator.items():
        # Compute across-all-samples major allele for one genomic site
        rcA, rcC, rcG, rcT, count_samples, scA, scC, scG, scT = site_info[:9]

        # Skip site with low prevalence for core sites and vice versa for rare sites
        prevalence = count_samples / total_samples_count
        if global_args.site_type == "common" and prevalence < global_args.site_prev:
            continue
        if global_args.site_type == "rare" and prevalence > global_args.site_prev:
            continue

        # Compute the pooled major allele based on the pooled-read-counts (abundance) or pooled-sample-counts (prevalence)
        if global_args.snp_pooled_method == "abundance":
            site_depth = sum((rcA, rcC, rcG, rcT))
            tuple_of_alleles = (('A', rcA), ('C', rcC), ('G', rcG), ('T', rcT))
        else:
            site_depth = count_samples
            tuple_of_alleles = (('A', scA), ('C', scC), ('G', scG), ('T', scT))

        major_allele, minor_allele, snp_type, number_alleles = call_alleles(tuple_of_alleles, site_depth, global_args.snp_maf)

        # Keep sites with desired snp_type
        if ('any' not in global_args.snp_type and snp_type not in global_args.snp_type):
            continue

        if number_alleles == 0:
            continue


        major_index = 'ACGT'.index(major_allele)
        minor_index = 'ACGT'.index(minor_allele)

        # Extract the read counts of previously computed across-samples major alleles
        sample_depths = [] # only accounts for reads matching either major or minor allele
        sample_mafs = [] # frequency of minor allele frequency
        for sample_index in range(9, len(site_info)):
            # For Each <Site, Sample> Pair
            rc_ACGT = [int(rc) for rc in site_info[sample_index].split(",")]

            sample_depth = rc_ACGT[major_index] if major_index == minor_index else rc_ACGT[major_index] + rc_ACGT[minor_index]
            maf_by_sample = -1.0 if sample_depth == 0 else (0.0 if major_index == minor_index else rc_ACGT[minor_index] / sample_depth)

            sample_depths.append(sample_depth)
            sample_mafs.append(maf_by_sample)

        # Site Annotation
        ref_id, ref_pos, ref_allele = site_id.rsplit("|", 2)
        ref_pos = int(ref_pos) # ref_pos is 1-based
        annots = ("IGR",) # short contigs may not carry any gene
        if ref_id in genes_boundary:
            annots = annotate_site(ref_pos, genes_boundary[ref_id], genes_feature[ref_id], genes_sequence) #<--

        locus_type = annots[0]
        gene_id = annots[1] if len(annots) > 1 else None
        site_type = annots[2] if len(annots) > 2 else None
        amino_acids = annots[3] if len(annots) > 2 else None

        pooled_snps_dict["info"][site_id] = {
            "site_id": site_id,
            "major_allele": major_allele,
            "minor_allele": minor_allele,
            "count_samples": count_samples,
            "snp_type": snp_type,
            "rcA": rcA, "rcC": rcC, "rcG": rcG, "rcT": rcT,
            "scA": scA, "scC": scC, "scG": scG, "scT":scT,
            "locus_type": locus_type,
            "gene_id": gene_id,
            "site_type": site_type,
            "amino_acids": amino_acids
        }
        pooled_snps_dict["freq"][site_id] = [site_id] + sample_mafs
        pooled_snps_dict["depth"][site_id] = [site_id] + sample_depths

    return pooled_snps_dict


def write_population_snps(pooled_snps_dict, species_id, chunk_id):

    global pool_of_samples

    if chunk_id == -2:
        # species
        global dict_of_species
        samples_names = dict_of_species[species_id].fetch_samples_names()
        snps_info_fp = pool_of_samples.get_target_layout("snps_info", species_id)
        snps_freq_fp = pool_of_samples.get_target_layout("snps_freq", species_id)
        snps_depth_fp = pool_of_samples.get_target_layout("snps_depth", species_id)
    else:
        # chunk
        snps_info_fp = pool_of_samples.get_target_layout("snps_info_by_chunk", species_id, chunk_id)
        snps_freq_fp = pool_of_samples.get_target_layout("snps_freq_by_chunk", species_id, chunk_id)
        snps_depth_fp = pool_of_samples.get_target_layout("snps_depth_by_chunk", species_id, chunk_id)

    with OutputStream(snps_info_fp) as stream:
        if chunk_id == -2:
            stream.write("\t".join(list(snps_info_schema.keys())) + "\n")
        for line in pooled_snps_dict["info"].values():
            stream.write("\t".join(map(format_data, line.values())) + "\n")

    with OutputStream(snps_freq_fp) as stream:
        if chunk_id == -2:
            stream.write("site_id\t" + "\t".join(samples_names) + "\n")
        for line in pooled_snps_dict["freq"].values():
            stream.write("\t".join(map(format_data, line)) + "\n")

    with OutputStream(snps_depth_fp) as stream:
        if chunk_id == -2:
            stream.write("site_id\t" + "\t".join(samples_names) + "\n")
        for line in pooled_snps_dict["depth"].values():
            stream.write("\t".join(map(str, line)) + "\n")


def collect_chunks(species_id):

    global global_args
    global dict_of_species
    global pool_of_samples

    sp = dict_of_species[species_id]
    number_of_chunks = sp.num_of_snps_chunks
    samples_names = dict_of_species[species_id].fetch_samples_names()

    loc_snps_info = [pool_of_samples.get_target_layout("snps_info_by_chunk", species_id, chunk_id) for chunk_id in range(0, number_of_chunks)]
    species_snps_info_fp = pool_of_samples.get_target_layout("snps_info", species_id)
    with OutputStream(species_snps_info_fp) as stream:
        stream.write("\t".join(list(snps_info_schema.keys())) + "\n")
    cat_files(loc_snps_info, species_snps_info_fp, 10)

    loc_snps_freq = [pool_of_samples.get_target_layout("snps_freq_by_chunk", species_id, chunk_id) for chunk_id in range(0, number_of_chunks)]
    species_snps_freq_fp = pool_of_samples.get_target_layout("snps_freq", species_id)
    with OutputStream(species_snps_freq_fp) as stream:
        stream.write("site_id\t" + "\t".join(samples_names) + "\n")
    cat_files(loc_snps_freq, species_snps_freq_fp, 10)

    loc_snps_depth = [pool_of_samples.get_target_layout("snps_depth_by_chunk", species_id, chunk_id) for chunk_id in range(0, number_of_chunks)]
    species_snps_depth = pool_of_samples.get_target_layout("snps_depth", species_id)
    with OutputStream(species_snps_depth) as stream:
        stream.write("site_id\t" + "\t".join(samples_names) + "\n")
    cat_files(loc_snps_depth, species_snps_depth, 10)

    if not global_args.debug:
        for s_file in loc_snps_info + loc_snps_freq + loc_snps_depth:
            command(f"rm -rf {s_file}", quiet=True)
    return True


def merge_snps(args):

    try:
        global global_args
        global_args = args

        global pool_of_samples
        global dict_of_species

        pool_of_samples = SamplePool(args.samples_list, args.midas_outdir, "snps")
        assert len(pool_of_samples.samples) > 0, f"No samples in the provided samples_list"

        dict_of_species = pool_of_samples.select_species("snps", args)
        species_ids_of_interest = [sp.id for sp in dict_of_species.values()]
        assert species_ids_of_interest, f"No (specified) species pass the genome_coverage filter across samples, please adjust the genome_coverage or species_list"
        species_count = len(species_ids_of_interest)
        tsprint(f"{species_count} species pass the filter")

        pool_of_samples.create_dirs(["outdir", "tempdir"], args.debug)
        pool_of_samples.create_species_subdirs(species_ids_of_interest, "outdir", args.debug, quiet=True)
        pool_of_samples.create_species_subdirs(species_ids_of_interest, "tempdir", args.debug, quiet=True)

        tsprint(f"MIDAS2::write_species_summary::start")
        pool_of_samples.write_summary_files(dict_of_species, "snps")
        tsprint(f"MIDAS2::write_species_summary::finish")

        # Download representative genomes for every species into midas_db
        num_cores_download = min(args.num_cores, len(species_ids_of_interest))
        midas_db = MIDAS_DB(os.path.abspath(args.midasdb_dir), args.midasdb_name, num_cores_download)

        # The unit of compute across-samples pop SNPs is: chunk_of_sites.
        tsprint(f"MIDAS2::design_chunks::start")
        arguments_list = design_chunks(species_ids_of_interest, midas_db)
        tsprint(f"MIDAS2::design_chunks::finish")

        tsprint(f"MIDAS2::multiprocessing_map::start")
        proc_flags = multiprocessing_map(process, arguments_list, args.num_cores)
        assert all(s == "worked" for s in proc_flags), f"Error: some chunks failed"
        tsprint(f"MIDAS2::multiprocessing_map::finish")

    except AssertionError as error:
        tsprint(f"Bugs in the codes, keep the outputs for debugging purpose.")
        # ideally exiting the program, logging the issue, and notifying the user.
        raise error
    except Exception as error:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error. Specify --debug flag to keep.")
            pool_of_samples.remove_dirs(["outdir", "tempdir"])
        raise error


@register_args
def main(args):
    tsprint(f"Population SNPs calling in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    merge_snps(args)
