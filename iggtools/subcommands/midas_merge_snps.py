
import json
from collections import defaultdict
from operator import itemgetter
import multiprocessing
from math import ceil
from bisect import bisect
import Bio.SeqIO

from iggtools.models.samplepool import SamplePool
from iggtools.common.utils import tsprint, num_physical_cores, command, InputStream, OutputStream, multiprocessing_map, multithreading_map, select_from_tsv
from iggtools.models.uhgg import MIDAS_IGGDB
from iggtools.params.schemas import snps_pileup_schema, snps_info_schema, format_data, genes_feature_schema
from iggtools.subcommands.midas_run_snps import cat_files, scan_contigs
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


def fetch_ref_codon(ref_pos, curr_gene, curr_seq):
    """ Fetch codon within gene for given site """
    # position of site in gene
    within_gene_position = ref_pos - curr_gene['start'] if curr_gene['strand'] == '+' else curr_gene['end'] - ref_pos
    # position of site in codon
    within_codon_position = within_gene_position % 3
    # gene sequence (oriented start to stop)
    ref_codon = curr_seq[within_codon_position-within_codon_position:within_codon_position-within_codon_position+3]
    return ref_codon, within_codon_position


def read_gene_sequence(fasta_file):
    """ Scan the genome file to get contig_id and contig_seq as ref_seq """
    contigs = {}
    with InputStream(fasta_file) as file:
        for rec in Bio.SeqIO.parse(file, 'fasta'):
            contigs[rec.id] = {
                "gene_id": rec.id,
                "gene_len": len(rec.seq),
                "gene_seq": str(rec.seq),
            }
    ## need to pay attention to reverse strand, probably by using scan_contigs()
    return contigs


def read_gene_features(features_file):
    """ Read TAB-delimited *.genes files from gene_annotations """
    features = defaultdict(dict)
    with InputStream(features_file) as stream:
        for r in select_from_tsv(stream, selected_columns=genes_feature_schema, result_structure=dict):
            cid = r['contig_id']
            r['contig_id'] = f"gnl|Prokka|{cid}" ## TODO: build gene feature should not have replaced the prefix in the first place
            features[r['contig_id']][r['gene_id']] = r
    return features


def check_feature_counts(features, gene_seqs):
    """ Check if the parsed gene feature file is consistent with the Prokka gene ffn file """
    counts = 0
    genes = []
    for _ in features.keys():
        counts += len(features[_])
        genes.extend(list(features[_].keys()))
    assert len(gene_seqs) == counts
    return sorted(genes) == sorted(list(gene_seqs.keys()))


def check_gene_sequences(contig_file, features, gene_seqs, species_id):
    """ Check if the prokka gene sequences is same with extracting from genome"""
    contig_seqs = scan_contigs(contig_file, species_id)
    flags = dict()
    for ref_id, c in contig_seqs.items():
        f = features[ref_id]
        for gid, gdict in f.items():
            flags[gid] = gene_seqs[gid]["gene_seq"] == get_gen_seq(c["contig_seq"], gdict['start'], gdict['end'], gdict['strand'])
    return all(list(flags.values()))


def binary_search_site(list_of_boundaries, ref_pos):
    """ Binary search the boundaries, if return odd than within-ranges otherwise between-ranges """
    flag = bisect(list_of_boundaries, ref_pos)
    if flag % 2: # even: intergenic
        return None
    index = int((flag + 1) / 2)
    ## Return the index of the ranges (1-based)
    return index


def generate_boundaries(features):
    """ Given list of gene ranges, generate the desired, half-open boundaries by binary search """
    gene_boundaries = defaultdict(dict)
    for contig_id in features.keys():
        feature_per_contig = features[contig_id]
        # Sort features by starting position
        feature_per_contig_sorted = dict(sorted(feature_per_contig.items(), key=feature_per_contig.get("start"), reverse=False))
        feature_ranges = {gf['gene_id']: (gf['start'], gf['end']) for gf in feature_per_contig_sorted.values()} ## +1 double check
        # TODO: double check non-overlapping ranges before flatten
        feature_ranges_flat = tuple(_ for rt in tuple(feature_ranges.values()) for _ in rt)
        # Convert ranges into half-open intervals.
        boundaries = tuple(gr + 1 if idx%2 == 1 else gr for idx, gr in enumerate(feature_ranges_flat))
        gene_boundaries[contig_id] = {"genes": list(feature_ranges.keys()), "boundaries": boundaries} #, "ranges": feature_ranges}
        ## REMOVE ranges later after debugging mode
    return gene_boundaries


def compute_degenracy(ref_codon, within_codon_position, strand):
    """ Compute degenracy """
    amino_acids = []
    for allele in ['A', 'C', 'G', 'T']: # + strand
        codon = index_replace(ref_codon, allele, within_codon_position, strand) # +/- strand
        amino_acid = translate(codon)
        amino_acids.append(amino_acid)
    unique_aa = set(amino_acids)
    degeneracy = 4 - len(unique_aa) + 1
    site_type = f"{degeneracy}D"
    amino_acids = ','.join(amino_acids)
    return site_type, amino_acids


def annotate_site(ref_id, ref_pos, curr_contig, curr_feature, gene_seqs):
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

    curr_seq = gene_seqs[curr_gene_id]["gene_seq"]
    assert len(curr_seq) % 3 == 0, f"gene must by divisible by 3 to id codons"
    ref_codon, within_codon_position = fetch_ref_codon(ref_pos, curr_gene, curr_seq)
    assert all(_ in ['A', 'T', 'C', 'G'] for _ in ref_codon), f"codon {ref_codon} for {ref_id}-{ref_pos} contain weird characters"

    site_type, amino_acids = compute_degenracy(ref_codon, within_codon_position, curr_gene['strand'])
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


def design_chunks(contigs_files, gene_features_files, gene_seqs_files, chunk_size):
    """ Chunks_of_continuous_genomic_sites and each chunk is indexed by (species_id, chunk_id) """

    tsprint("CZ::design_chunks::start species_counts: %s" % (len(contigs_files)))

    # dict of semaphore and acquire before
    global semaphore_for_species
    global species_sliced_pileup_path
    global species_samples_dict

    global pool_of_samples
    global dict_of_species

    semaphore_for_species = dict()
    species_sliced_pileup_path = defaultdict(dict)
    species_samples_dict = defaultdict(dict)

    argument_list = []
    for species in dict_of_species.values():
        species_id = species.id
        contigs = scan_contigs(contigs_files[species_id], species_id)

        samples_depth = species.samples_depth
        samples_snps_pileup = [sample.get_target_layout("snps_pileup", species_id) for sample in list(species.samples)]

        species_samples_dict["samples_depth"][species_id] = samples_depth
        species_samples_dict["samples_snps_pileup"][species_id] = samples_snps_pileup

        total_samples_count = len(species.samples)

        chunk_id = 0
        for contig_id, contig in contigs.items():
            contig_length = contig["contig_len"]

            # pileup is 1-based index
            if contig_length <= chunk_size:
                my_args = (species_id, chunk_id, contig_id, 1, contig_length, total_samples_count, gene_features_files[species_id], gene_seqs_files[species_id]) ##cz
                argument_list.append(my_args)

                snps_info_fp = pool_of_samples.get_target_layout("snps_info_by_chunk", species_id, chunk_id)
                snps_freq_fp = pool_of_samples.get_target_layout("snps_freq_by_chunk", species_id, chunk_id)
                snps_depth_fp = pool_of_samples.get_target_layout("snps_depth_by_chunk", species_id, chunk_id)
                species_sliced_pileup_path[species_id][chunk_id] = (snps_info_fp, snps_freq_fp, snps_depth_fp)
                chunk_id += 1
            else:
                number_of_chunks = ceil(contig_length/chunk_size) - 1
                for ni, ci in enumerate(range(0, contig_length, chunk_size)):
                    if ni == number_of_chunks:
                        my_args = (species_id, chunk_id, contig_id, ci+1, contig_length, total_samples_count, gene_features_files[species_id], gene_seqs_files[species_id]) ##cz
                    else:
                        my_args = (species_id, chunk_id, contig_id, ci+1, ci+chunk_size, total_samples_count, gene_features_files[species_id], gene_seqs_files[species_id]) ##cz
                    argument_list.append(my_args)

                    snps_info_fp = pool_of_samples.get_target_layout("snps_info_by_chunk", species_id, chunk_id)
                    snps_freq_fp = pool_of_samples.get_target_layout("snps_freq_by_chunk", species_id, chunk_id)
                    snps_depth_fp = pool_of_samples.get_target_layout("snps_depth_by_chunk", species_id, chunk_id)
                    species_sliced_pileup_path[species_id][chunk_id] = (snps_info_fp, snps_freq_fp, snps_depth_fp)
                    chunk_id += 1
        tsprint(f"  CZ::design_chunks::{species_id}::finish for loop with {chunk_id} chunks")

        # Submit the merge jobs
        argument_list.append((species_id, -1))
        snps_info_fp = pool_of_samples.get_target_layout("snps_info", species_id)
        snps_freq_fp = pool_of_samples.get_target_layout("snps_freq", species_id)
        snps_depth_fp = pool_of_samples.get_target_layout("snps_depth", species_id)
        species_sliced_pileup_path[-1][species_id] = (snps_info_fp, snps_freq_fp, snps_depth_fp)

        # Create a semaphore with number_of_chunks for current species
        semaphore_for_species[species_id] = multiprocessing.Semaphore(chunk_id)
        for _ in range(chunk_id):
            semaphore_for_species[species_id].acquire()

    tsprint("CZ::design_chunks::finish species_counts: %s" % (len(contigs_files)))
    return argument_list


def merge_chunks_by_species(species_id):

    tsprint(f"    CZ::merge_chunks_by_species::{species_id}::start")

    global global_args
    global species_sliced_pileup_path
    global dict_of_species

    snps_info_fp, snps_freq_fp, snps_depth_fp = species_sliced_pileup_path[-1][species_id]
    sliced_pileup_files = list(species_sliced_pileup_path[species_id].values())
    snps_info_files = [fi[0] for fi in sliced_pileup_files]
    snps_freq_files = [fi[1] for fi in sliced_pileup_files]
    snps_depth_files = [fi[2] for fi in sliced_pileup_files]
    samples_names = dict_of_species[species_id].fetch_samples_names()

    # Add header for the merged-chunks
    with OutputStream(snps_info_fp) as out_info:
        out_info.write("\t".join(list(snps_info_schema.keys())) + "\n")
    cat_files(snps_info_files, snps_info_fp, 20)

    with OutputStream(snps_freq_fp) as out_freq:
        out_freq.write("site_id\t" + "\t".join(samples_names) + "\n")
    cat_files(snps_freq_files, snps_freq_fp, 10)

    with OutputStream(snps_depth_fp) as out_depth:
        out_depth.write("site_id\t" + "\t".join(samples_names) + "\n")
    cat_files(snps_depth_files, snps_depth_fp, 10)

    if not global_args.debug:
        for s_file in snps_info_files + snps_freq_files + snps_depth_files:
            command(f"rm -rf {s_file}", quiet=True)

    tsprint(f"    CZ::merge_chunks_by_species::{species_id}::finish")
    return True


def accumulate(accumulator, proc_args):
    """ Accumulate read_counts and sample_counts for a chunk of sites for one sample,
    at the same time remember <site, sample>'s A, C, G, T read counts."""

    contig_id, contig_start, contig_end, sample_index, snps_pileup_path, total_samples_count, genome_coverage = proc_args
    tsprint(f"    CZ::accumulate::{contig_id}-{contig_start}-{sample_index}::start")

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

    tsprint(f"    CZ::accumulate::{contig_id}-{contig_start}-{sample_index}::finish")


def compute_and_write_pooled_snps(accumulator, total_samples_count, species_id, chunk_id, gene_feature_file, gene_seq_file):
    """ For each site, compute the pooled-major-alleles, site_depth, and vector of sample_depths and sample_minor_allele_freq"""

    tsprint(f"    CZ::compute_and_write_pooled_snps::{species_id}-{chunk_id}::start")

    global global_args
    args = global_args

    global species_sliced_pileup_path
    snps_info_fp, snps_freq_fp, snps_depth_fp = species_sliced_pileup_path[species_id][chunk_id]

    # TODO: how to pass on these two dictionaries, instead in parse it over and over again in each chunk ...?
    features_by_contig = read_gene_features(gene_feature_file)
    gene_boundaries = generate_boundaries(features_by_contig)
    gene_seqs = read_gene_sequence(gene_seq_file)

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
                tuple_of_alleles = (('A', rcA), ('C', rcC), ('G', rcG), ('T', rcT))
            else:
                site_depth = count_samples
                tuple_of_alleles = (('A', scA), ('C', scC), ('G', scG), ('T', scT))

            major_allele, minor_allele, snp_type = call_alleles(tuple_of_alleles, site_depth, args.snp_maf)
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

                sample_depth = rc_ACGT[major_index] if major_index == minor_index else rc_ACGT[major_index] + rc_ACGT[minor_index]
                maf_by_sample = -1.0 if sample_depth == 0 else (0.0 if major_index == minor_index else rc_ACGT[minor_index] / sample_depth)

                sample_depths.append(sample_depth)
                sample_mafs.append(maf_by_sample)

            # Annotate one site
            ref_id, ref_pos, ref_allele = site_id.rsplit("|", 2)
            ref_pos = int(ref_pos) # ref_pos is 1-based
            if ref_id not in gene_boundaries:
                annots = ("IGR",) ## short contigs may not carry any gene
            else:
                curr_contig = gene_boundaries[ref_id]
                curr_feature = features_by_contig[ref_id]
                tsprint(f"    CZ::annotate_site::{ref_id}-{ref_pos}::start")
                annots = annotate_site(ref_id, ref_pos, curr_contig, curr_feature, gene_seqs)
                tsprint(f"    CZ::annotate_site::{ref_id}-{ref_pos}::finish {annots[0]}") # running time depends on the locus_type

            locus_type = annots[0]
            gene_id = annots[1] if len(annots) > 1 else None
            site_type = annots[2] if len(annots) > 2 else None
            amino_acids = annots[3] if len(annots) > 2 else None

            # Write
            out_info.write(f"{site_id}\t{major_allele}\t{minor_allele}\t{count_samples}\t{snp_type}\t{rcA}\t{rcC}\t{rcG}\t{rcT}\t{scA}\t{scC}\t{scG}\t{scT}\t{locus_type}\t{gene_id}\t{site_type}\t{amino_acids}\n")
            out_freq.write(f"{site_id}\t" + "\t".join(map(format_data, sample_mafs)) + "\n")
            out_depth.write(f"{site_id}\t" + "\t".join(map(str, sample_depths)) + "\n")

    tsprint(f"    CZ::compute_and_write_pooled_snps::{species_id}-{chunk_id}::finish")
    return True


def pool_one_chunk_across_samples(packed_args):
    """ For genome sites from one chunk, scan across all the sample, compute pooled SNPs and write to file """

    global semaphore_for_species
    global species_sliced_pileup_path
    global species_samples_dict

    species_id, chunk_id, contig_id, contig_start, contig_end, total_samples_count, gene_feature_file, gene_seq_file = packed_args
    tsprint(f"    CZ::pool_one_chunk_across_samples::{species_id}-{chunk_id}::start")

    list_of_snps_pileup_path = species_samples_dict["samples_snps_pileup"][species_id]
    list_of_sample_depths = species_samples_dict["samples_depth"][species_id]

    try:
        accumulator = dict()
        for sample_index in range(total_samples_count):
            proc_args = (contig_id, contig_start, contig_end, sample_index, list_of_snps_pileup_path[sample_index], total_samples_count, list_of_sample_depths[sample_index])
            accumulate(accumulator, proc_args)
        compute_and_write_pooled_snps(accumulator, total_samples_count, species_id, chunk_id, gene_feature_file, gene_seq_file)
        tsprint(f"    CZ::pool_one_chunk_across_samples::{species_id}-{chunk_id}::finish")
    finally:
        semaphore_for_species[species_id].release() # no deadlock


def process_chunk_of_sites(packed_args):

    global semaphore_for_species
    global species_sliced_pileup_path
    global species_samples_dict

    if packed_args[1] == -1:
        # Merge chunks_of_sites' pileup results per species
        species_id = packed_args[0]
        number_of_chunks = len(species_sliced_pileup_path[species_id])

        tsprint(f"  CZ::process_chunk_of_sites::{species_id}::wait merge_chunks_by_species")
        for _ in range(number_of_chunks):
            semaphore_for_species[species_id].acquire()
        tsprint(f"  CZ::process_chunk_of_sites::{species_id}::start merge_chunks_by_species")
        merge_chunks_by_species(species_id)
        tsprint(f"  CZ::process_chunk_of_sites::{species_id}::finish merge_chunks_by_species")
        return "worked"

    species_id, chunk_id = packed_args[:2]
    tsprint(f"  CZ::process_chunk_of_sites::{species_id}-{chunk_id}::start pool_one_chunk_across_samples")
    pool_one_chunk_across_samples(packed_args)
    tsprint(f"  CZ::process_chunk_of_sites::{species_id}-{chunk_id}::finish pool_one_chunk_across_samples")

    return "worked"


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
        midas_iggdb = MIDAS_IGGDB(args.midas_iggdb if args.midas_iggdb else pool_of_samples.get_target_layout("midas_iggdb_dir"), args.num_cores)
        contigs_files = midas_iggdb.fetch_files("prokka_genome", species_ids_of_interest) #contigs
        gene_features_files = midas_iggdb.fetch_files("gene_feature", species_ids_of_interest)
        gene_seqs_files = midas_iggdb.fetch_files("gene_seq", species_ids_of_interest)

        def check_annotation_setup(species_id):
            features_file = gene_features_files[species_id]
            gene_seq_file = gene_seqs_files[species_id]
            genome_file = contigs_files[species_id]

            features = read_gene_features(features_file)
            gene_seqs = read_gene_sequence(gene_seq_file)

            assert check_feature_counts(features, gene_seqs), f"Gene feature counts disagree with Prokka gene ffn file for species {species_id}"
            assert check_gene_sequences(genome_file, features, gene_seqs, species_id), f"Prokka gene sequences disagree with gene ranges computation for species {species_id}"

            return True

        assert all(multithreading_map(check_annotation_setup, species_ids_of_interest, num_threads=10))

        # Compute pooled SNPs by the unit of chunks_of_sites
        argument_list = design_chunks(contigs_files, gene_features_files, gene_seqs_files, args.chunk_size)
        proc_flags = multiprocessing_map(process_chunk_of_sites, argument_list, args.num_cores)
        assert all(s == "worked" for s in proc_flags)

    except Exception as error:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error. Specify --debug flag to keep.")
            pool_of_samples.remove_dirs(["outdir", "tempdir"])
        raise error


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_merge_snps(args)
