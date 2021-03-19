#!/usr/bin/env python3
import os
import json
from collections import defaultdict
from operator import itemgetter
import multiprocessing
from bisect import bisect
import time

from iggtools.models.samplepool import SamplePool
from iggtools.common.utils import tsprint, num_physical_cores, command, InputStream, OutputStream, multiprocessing_map, multithreading_map, select_from_tsv, cat_files
from iggtools.models.uhgg import MIDAS_IGGDB
from iggtools.params.schemas import snps_pileup_schema, snps_info_schema, format_data
from iggtools.models.species import Species
from iggtools.common.argparser import add_subcommand


DEFAULT_SAMPLE_COUNTS = 2
DEFAULT_GENOME_DEPTH = 5.0
DEFAULT_GENOME_COVERAGE = 0.4
DEFAULT_CHUNK_SIZE = 50000

DEFAULT_SITE_DEPTH = 2
DEFAULT_SITE_RATIO = 5.0

DEFAULT_SITE_PREV = 0.80
DEFAULT_SITE_TYPE = "common"

DEFAULT_SNP_POOLED_METHOD = "prevalence"
DEFAULT_SNP_MAF = 0.05
DEFAULT_SNP_TYPE = "mono, bi"


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
    subparser.add_argument('--chunk_size',
                           dest='chunk_size',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_CHUNK_SIZE,
                           help=f"Number of genomic sites for the temporary chunk file  ({DEFAULT_CHUNK_SIZE})")

    subparser.add_argument('--midas_iggdb',
                           dest='midas_iggdb',
                           type=str,
                           metavar="CHAR",
                           help=f"local MIDAS DB which mirrors the s3 IGG db")
    subparser.add_argument('--num_cores',
                           dest='num_cores',
                           type=int,
                           metavar="INT",
                           default=num_physical_cores,
                           help=f"Number of physical cores to use ({num_physical_cores})")

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

    return main_func


def acgt_string(A, C, G, T):
    return ','.join(map(str, (A, C, G, T)))


def translate(codon):
    """ Translate individual codon """
    codontable = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    return codontable[str(codon)]


def complement(base):
    """ Complement nucleotide """
    d = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    if base in d:
        return d[base]
    return base


def rev_comp(seq):
    """ Reverse complement sequence """
    return ''.join([complement(base) for base in list(seq[::-1])])


def get_gen_seq(genome_seq, start, end, strand):
    seq = genome_seq[start-1 : end]
    if strand == "-":
        return rev_comp(seq)
    return seq


def index_replace(codon, allele, pos, strand):
    """ Replace character at index i in string x with y"""
    bases = list(codon)
    bases[pos] = allele if strand == '+' else complement(allele)
    return ''.join(bases)


def binary_search_site(list_of_boundaries, ref_pos):
    """ Binary search the boundaries, if return odd than within-ranges otherwise between-ranges """
    flag = bisect(list_of_boundaries, ref_pos)
    if flag % 2 == 0: # even: intergenic
        return None
    index = int((flag + 1) / 2)
    # Return the index of the ranges (1-based)
    return index


def compute_degenracy(ref_codon, within_codon_pos, strand):
    """ Compute degenracy """
    amino_acids = []
    for allele in ['A', 'C', 'G', 'T']: # + strand
        codon = index_replace(ref_codon, allele, within_codon_pos, strand) # +/- strand
        amino_acid = translate(codon)
        amino_acids.append(amino_acid)
    unique_aa = set(amino_acids)
    degeneracy = 4 - len(unique_aa) + 1
    site_type = f"{degeneracy}D"
    amino_acids = ','.join(amino_acids)
    return site_type, amino_acids


def fetch_ref_codon(ref_pos, curr_gene, curr_seq):
    """ Fetch codon within gene for given site """
    # position of site in gene
    within_gene_pos = ref_pos - curr_gene['start'] if curr_gene['strand'] == '+' else curr_gene['end'] - ref_pos
    # position of site in codon
    within_codon_pos = within_gene_pos % 3
    # gene sequence (oriented start to stop)
    ref_codon = curr_seq[within_gene_pos-within_codon_pos:within_gene_pos-within_codon_pos+3]
    return ref_codon, within_codon_pos


def annotate_site(ref_id, ref_pos, curr_contig, curr_feature, genes_sequence):
    """ Annotate one genomic site, search against all genes for given species """
    # Binary search the range of the given genomic site position
    index = binary_search_site(curr_contig["boundaries"], ref_pos)
    if index is None:
        locus_type = "IGR" # even: intergenic
        return (locus_type,)

    curr_gene_id = curr_contig["genes"][index-1]
    curr_gene = curr_feature[curr_gene_id]
    locus_type = curr_gene["gene_type"]

    if locus_type != "CDS":
        return locus_type, curr_gene_id

    curr_seq = genes_sequence[curr_gene_id]["seq"]
    assert len(curr_seq) % 3 == 0, f"gene must by divisible by 3 to id codons"

    ref_codon, within_codon_pos = fetch_ref_codon(ref_pos, curr_gene, curr_seq)
    assert all(_ in ['A', 'T', 'C', 'G'] for _ in ref_codon), f"codon {ref_codon} for {ref_id}-{ref_pos} contain weird characters"

    site_type, amino_acids = compute_degenracy(ref_codon, within_codon_pos, curr_gene['strand'])
    return locus_type, curr_gene_id, site_type, amino_acids


def call_alleles(tuple_of_alleles, site_depth, snp_maf):
    """ Compute the pooled allele frequencies and call SNPs """

    # Only when you have seen all the revelant samples, you can call SNPs
    # keep alleles passing the min allele frequency
    alleles_above_cutoff = tuple(al for al in tuple_of_alleles if al[1] / site_depth >= snp_maf)

    # classify SNPs type
    number_alleles = len(alleles_above_cutoff)
    snp_type = ["mono", "bi", "tri", "quad"][number_alleles - 1]

    # In the event of a tie -- biallelic site with 50/50 freq split -- the allele declared major is
    # the one that comes later in the "ACGT" lexicographic order.
    alleles_above_cutoff = sorted(alleles_above_cutoff, key=itemgetter(1), reverse=True)[:2]
    major_allele = alleles_above_cutoff[0][0]
    minor_allele = alleles_above_cutoff[-1][0] # for fixed sites, same as major allele

    return (major_allele, minor_allele, snp_type)


def design_chunks(midas_iggdb, chunk_size):
    global pool_of_samples
    global dict_of_species

    global semaphore_for_species
    semaphore_for_species = dict()

    arguments_list = []
    for sp in dict_of_species.values():
        species_id = sp.id

        # The structure of the chunks depends on the representative genome sequences
        assert sp.design_snps_chunks(midas_iggdb, chunk_size)
        num_of_chunks = sp.num_of_sites_chunks

        for chunk_id in range(0, num_of_chunks):
            arguments_list.append((species_id, chunk_id))
        arguments_list.append((species_id, -1))

        semaphore_for_species[species_id] = multiprocessing.Semaphore(num_of_chunks)
        for _ in range(num_of_chunks):
            semaphore_for_species[species_id].acquire()

    return arguments_list


def prepare_annotation_per_species(args):
    sp, genes_feature_file, genes_seq_file = args
    sp.prepare_annotation(genes_feature_file, genes_seq_file)


def prepare_site_annotation(midas_iggdb, num_cores):
    global dict_of_species
    args_list = []
    for sp in dict_of_species.values():
        species_id = sp.id
        genes_feature_file = midas_iggdb.fetch_files("gene_feature", [species_id])[species_id]
        genes_seq_file = midas_iggdb.fetch_files("gene_seq", [species_id])[species_id]
        args_list.append((sp, genes_feature_file, genes_seq_file))
    multithreading_map(prepare_annotation_per_species, args_list, num_cores)


def process_one_chunk_of_sites(packed_args):

    species_id, chunk_id = packed_args

    if chunk_id == -1:
        global semaphore_for_species
        global dict_of_species
        sp = dict_of_species[species_id]

        tsprint(f"  CZ::process_one_chunk_of_sites::{species_id}--1::wait merge_all_chunks_per_species")
        for _ in range(sp.num_of_sites_chunks):
            semaphore_for_species[species_id].acquire()
        tsprint(f"  CZ::process_one_chunk_of_sites::{species_id}--1::start merge_all_chunks_per_species")
        merge_all_chunks_per_species(species_id)
        tsprint(f"  CZ::process_one_chunk_of_sites::{species_id}--1::finish merge_all_chunks_per_species")
        return "worked"

    tsprint(f"  CZ::process_one_chunk_of_sites::{species_id}-{chunk_id}::start pool_across_samples_per_chunk")
    pool_across_samples_per_chunk(species_id, chunk_id)
    tsprint(f"  CZ::process_one_chunk_of_sites::{species_id}-{chunk_id}::finish pool_across_samples_per_chunk")

    return "worked"


def pool_across_samples_per_chunk(species_id, chunk_id):
    """ For genome sites from one chunk, scan across all the sample, compute pooled SNPs and write to file """

    global semaphore_for_species
    global pool_of_samples
    global dict_of_species

    try:
        sp = dict_of_species[species_id]
        tsprint(f"    CZ::pool_across_samples_per_chunk::{species_id}-{chunk_id}::start accumulate_samples_per_unit")
        flags = []
        for pargs in sp.chunks_of_sites[chunk_id]:
            flags.append(accumulate_samples_per_unit(pargs))
        assert all(flags)
        tsprint(f"    CZ::pool_across_samples_per_chunk::{species_id}-{chunk_id}::finish accumulate_samples_per_unit")
    finally:
        semaphore_for_species[species_id].release() # no deadlock


def accumulate_samples_per_unit(packed_args):
    global pool_of_samples
    global dict_of_species

    species_id, chunk_id, contig_id, contig_start, contig_end = packed_args[:5]

    sp = dict_of_species[species_id]
    total_samples_count = sp.samples_count
    list_of_samples_depth = sp.list_of_samples_depth

    # Accumulate sites sample by sample
    tsprint(f"    CZ::accumulate_samples_per_unit::{species_id}-{chunk_id}::start accumulate_sample_by_sample")
    accumulator = dict()
    for sample_index, sample in enumerate(sp.list_of_samples):
        species_pileup_path = sample.get_target_layout("snps_pileup", species_id)
        chunk_pileup_path = sample.get_target_layout("chunk_pileup", species_id, chunk_id) # USE headerless_chunk_pileup_file if exits
        snps_pileup_path = chunk_pileup_path if os.path.exists(chunk_pileup_path) else species_pileup_path
        has_header = not os.path.exists(chunk_coverage_path)

        # Pileup is 1-based index, close left close right
        proc_args = (contig_id, contig_start+1, contig_end, sample_index, snps_pileup_path, total_samples_count, list_of_samples_depth[sample_index])
        accumulate(accumulator, proc_args)
    tsprint(f"    CZ::accumulate_samples_per_unit::{species_id}-{chunk_id}::finish accumulate_one_samples")

    # Compute across-samples SNPs and write to chunk file
    tsprint(f"    CZ::accumulate_samples_per_unit::{species_id}-{chunk_id}::start compute_and_write_pooled_snps_per_unit")
    flag = compute_and_write_pooled_snps_per_unit(accumulator, species_id, chunk_id)
    tsprint(f"    CZ::accumulate_samples_per_unit::{species_id}-{chunk_id}::finish compute_and_write_pooled_snps_per_unit")

    return flag


def accumulate(accumulator, proc_args):
    """ Accumulate read_counts and sample_counts for a chunk of sites for one sample,
    at the same time remember <site, sample>'s A, C, G, T read counts."""

    global global_args

    contig_id, contig_start, contig_end, sample_index, snps_pileup_path, total_samples_count, genome_coverage = proc_args

    # Output column indices
    c_A, c_C, c_G, c_T, c_count_samples, c_scA, c_scC, c_scG, c_scT = range(9)

    awk_command = f"awk \'$1 == \"{contig_id}\" && $2 >= {contig_start} && $2 <= {contig_end}\'"
    with InputStream(snps_pileup_path, awk_command) as stream:
        for row in select_from_tsv(stream, schema=snps_pileup_schema, result_structure=dict):
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
            assert acc[9 + sample_index] == '0,0,0,0' and acgt_str != '0,0,0,0', f"accumulate error::{site_id}:{acc}:{sample_index}"
            acc[9 + sample_index] = acgt_str


def compute_and_write_pooled_snps_per_unit(accumulator, species_id, chunk_id):
    """ For each site, compute the pooled-major-alleles, site_depth, and vector of sample_depths and sample_minor_allele_freq"""
    global global_args
    global dict_of_species

    sp = dict_of_species[species_id]
    total_samples_count = sp.samples_count

    pooled_snps_info_dict = {}
    pooled_snps_freq_list = []
    pooled_snps_depth_list = []

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

        major_allele, minor_allele, snp_type = call_alleles(tuple_of_alleles, site_depth, global_args.snp_maf)
        major_index = 'ACGT'.index(major_allele)
        minor_index = 'ACGT'.index(minor_allele)

        # Keep sites with desired snp_type
        if ('any' not in global_args.snp_type and snp_type not in global_args.snp_type):
            continue

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
        annots = ("IGR",)
        if ref_id not in sp.genes_boundary:
            annots = ("IGR",) # short contigs may not carry any gene
        else:
            annots = annotate_site(ref_id, ref_pos, sp.genes_boundary[ref_id], sp.genes_feature[ref_id], sp.genes_sequence)

        locus_type = annots[0]
        gene_id = annots[1] if len(annots) > 1 else None
        site_type = annots[2] if len(annots) > 2 else None
        amino_acids = annots[3] if len(annots) > 2 else None

        pooled_snps_info_dict[site_id] = {"site_id": site_id,
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
        pooled_snps_freq_list.append([site_id] + sample_mafs)
        pooled_snps_depth_list.append([site_id] + sample_depths)

    # Write to file
    snps_info_fp = pool_of_samples.get_target_layout("snps_info_by_chunk", species_id, chunk_id)
    snps_freq_fp = pool_of_samples.get_target_layout("snps_freq_by_chunk", species_id, chunk_id)
    snps_depth_fp = pool_of_samples.get_target_layout("snps_depth_by_chunk", species_id, chunk_id)

    with open(snps_info_fp, "a") as out_info:
        for r in pooled_snps_info_dict.values():
            out_info.write("\t".join(map(format_data, r.values())) + "\n")

    with open(snps_freq_fp, "a") as out_freq:
        for line in pooled_snps_freq_list:
            out_freq.write("\t".join(map(format_data, line)) + "\n")

    with open(snps_depth_fp, "a") as out_depth:
        for line in pooled_snps_depth_list:
            out_depth.write("\t".join(map(str, line)) + "\n")

    return True


def merge_all_chunks_per_species(species_id):

    global global_args
    global dict_of_species
    global pool_of_samples

    sp = dict_of_species[species_id]
    number_of_chunks = sp.num_of_sites_chunks
    samples_names = dict_of_species[species_id].fetch_samples_names()

    list_of_chunks_snps_info = [pool_of_samples.get_target_layout("snps_info_by_chunk", species_id, chunk_id) for chunk_id in range(0, number_of_chunks)]
    species_snps_info_fp = pool_of_samples.get_target_layout("snps_info", species_id)

    # Add header
    with OutputStream(species_snps_info_fp) as stream:
        stream.write("\t".join(list(snps_info_schema.keys())) + "\n")
    cat_files(list_of_chunks_snps_info, species_snps_info_fp, 10)

    list_of_chunks_snps_freq = [pool_of_samples.get_target_layout("snps_freq_by_chunk", species_id, chunk_id) for chunk_id in range(0, number_of_chunks)]
    species_snps_freq_fp = pool_of_samples.get_target_layout("snps_freq", species_id)
    with OutputStream(species_snps_freq_fp) as stream:
        stream.write("site_id\t" + "\t".join(samples_names) + "\n")
    cat_files(list_of_chunks_snps_freq, species_snps_freq_fp, 10)

    list_of_chunks_snps_depth = [pool_of_samples.get_target_layout("snps_depth_by_chunk", species_id, chunk_id) for chunk_id in range(0, number_of_chunks)]
    species_snps_depth = pool_of_samples.get_target_layout("snps_depth", species_id)
    with OutputStream(species_snps_depth) as stream:
        stream.write("site_id\t" + "\t".join(samples_names) + "\n")
    cat_files(list_of_chunks_snps_depth, species_snps_depth, 10)

    #if not global_args.debug:
    #    for s_file in list_of_chunks_snps_info + list_of_chunks_snps_freq + list_of_chunks_snps_depth:
    #        command(f"rm -rf {s_file}", quiet=True)

    return True


def midas_merge_snps(args):

    try:
        global global_args
        global_args = args

        global pool_of_samples
        global dict_of_species

        pool_of_samples = SamplePool(args.samples_list, args.midas_outdir, "snps")
        dict_of_species = pool_of_samples.select_species("snps", args)

        species_ids_of_interest = [sp.id for sp in dict_of_species.values()]
        assert species_ids_of_interest, f"No (specified) species pass the genome_coverage filter across samples, please adjust the genome_coverage or species_list"
        tsprint(species_ids_of_interest)


        pool_of_samples.create_dirs(["outdir", "tempdir"], args.debug)
        pool_of_samples.create_species_subdirs(species_ids_of_interest, "outdir", args.debug)
        pool_of_samples.create_species_subdirs(species_ids_of_interest, "tempdir", args.debug)
        pool_of_samples.write_summary_files(dict_of_species, "snps")


        # Download representative genomes for every species into midas_iggdb
        num_cores = min(args.num_cores, len(species_ids_of_interest))
        midas_iggdb = MIDAS_IGGDB(args.midas_iggdb if args.midas_iggdb else pool_of_samples.get_target_layout("midas_iggdb_dir"), num_cores)

        # The unit of compute across-samples pooled SNPs is: chunk_of_sites.
        tsprint(f"CZ::design_chunks::start")
        arguments_list = design_chunks(midas_iggdb, args.chunk_size)
        tsprint(f"CZ::design_chunks::finish")


        tsprint(f"CZ::prepare_site_annotation::start")
        prepare_site_annotation(midas_iggdb, num_cores)
        tsprint(f"CZ::prepare_site_annotation::finish")


        tsprint(f"CZ::multiprocessing_map::start")
        proc_flags = multiprocessing_map(process_one_chunk_of_sites, arguments_list, args.num_cores)
        tsprint(f"CZ::multiprocessing_map::finish")

        assert all(s == "worked" for s in proc_flags)

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
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_merge_snps(args)
