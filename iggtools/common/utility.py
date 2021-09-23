#!/usr/bin/env python3
from bisect import bisect


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


def annotate_site(ref_pos, curr_contig, curr_feature, genes_sequence):
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

    # Codon contains weird characters.
    if not all(_ in ['A', 'T', 'C', 'G'] for _ in ref_codon):
        return locus_type, curr_gene_id

    site_type, amino_acids = compute_degenracy(ref_codon, within_codon_pos, curr_gene['strand'])
    return locus_type, curr_gene_id, site_type, amino_acids
