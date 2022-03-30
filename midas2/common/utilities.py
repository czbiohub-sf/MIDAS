#!/usr/bin/env python3
import os
from bisect import bisect
from collections import defaultdict
import Bio.SeqIO

from midas2.common.utils import InputStream, retry, select_from_tsv, tsprint
from midas2.params.schemas import genes_feature_schema, PAN_GENE_INFO_SCHEMA, MARKER_INFO_SCHEMA, PAN_GENE_LENGTH_SCHEMA, CLUSTER_INFO_SCHEMA


def decode_species_arg(args, species):
    selected_species = set()
    try:  # pylint: disable=too-many-nested-blocks
        if args.species.upper() == "ALL":
            selected_species = set(species)
        else:
            for s in args.species.split(","):
                if ":" not in s:
                    assert str(int(s)) == s, f"Species id is not an integer: {s}"
                    selected_species.add(s)
                else:
                    i, n = s.split(":")
                    i = int(i)
                    n = int(n)
                    assert 0 <= i < n, f"Species class and modulus make no sense: {i}, {n}"
                    for sid in species:
                        if int(sid) % n == i:
                            selected_species.add(sid)
    except:
        tsprint(f"ERROR:  Species argument is not a list of species ids or slices: {s}")
        raise
    return sorted(selected_species)


def decode_genomes_arg(args, genomes):
    selected_genomes = set()
    try:  # pylint: disable=too-many-nested-blocks
        if args.genomes.upper() == "ALL":
            selected_genomes = set(genomes)
        else:
            for g in args.genomes.split(","):
                if ":" not in g:
                    selected_genomes.add(g)
                else:
                    i, n = g.split(":")
                    i = int(i)
                    n = int(n)
                    assert 0 <= i < n, f"Genome class and modulus make no sense: {i}, {n}"
                    for gid in genomes:
                        gid_int = int(gid.replace("GUT_GENOME", ""))
                        if gid_int % n == i:
                            selected_genomes.add(gid)
    except:
        tsprint(f"ERROR:  Genomes argument is not a list of genome ids or slices: {g}")
        raise
    return sorted(selected_genomes)


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



@retry
def scan_fasta(fasta_file):
    """ Scan FASTA FILE to get seq and len """
    seqs = {}
    with InputStream(fasta_file) as file:
        for rec in Bio.SeqIO.parse(file, 'fasta'):
            seqs[rec.id] = {
                "id": rec.id,
                "length": len(rec.seq),
                "seq": str(rec.seq),
            }
    return seqs


def update_id(cid):
    cid = cid.replace("gnl|Prokka|", "") #<-------
    cid = cid.replace("UHGGGC", "GC")
    return cid


@retry
def scan_mapfile(mapfile):
    """ Extract <marker_id, gene_id> pairs from the marker mapfile """
    dict_of_markers = dict()
    with InputStream(mapfile) as stream:
        for gene_id, marker_id in select_from_tsv(stream, ["gene_id", "marker_id"], schema=MARKER_INFO_SCHEMA):
            dict_of_markers[gene_id] = marker_id
    return dict_of_markers


@retry
def scan_gene_info(gene_info_fp):
    """ Extract <centroid_99, ...> from the gene_info file """
    dict_of_centroids = dict()
    with InputStream(gene_info_fp) as stream:
        for r in select_from_tsv(stream, selected_columns=list(PAN_GENE_INFO_SCHEMA.keys())[1:], result_structure=dict):
            if r["centroid_99"] not in dict_of_centroids:
                dict_of_centroids[r["centroid_99"]] = r
    return dict_of_centroids


@retry
def scan_gene_length(gene_length_file):
    gene_length_dict = dict()
    with InputStream(gene_length_file) as stream:
        for r in select_from_tsv(stream, schema=PAN_GENE_LENGTH_SCHEMA, result_structure=dict):
            gene_length_dict[r["gene_id"]] = r["gene_length"]
    return gene_length_dict


@retry
def scan_cluster_info(info_file):
    """ Read in cluster_info """
    cluster_info = dict()
    with InputStream(info_file) as stream:
        for r in select_from_tsv(stream, selected_columns=CLUSTER_INFO_SCHEMA, result_structure=dict):
            cluster_info[r["centroid_99"]] = r
    return cluster_info


@retry
def scan_genes(annotated_genes):
    """" Lookup of seq_id to sequence for PATRIC genes """
    gene_seqs = {}
    with InputStream(annotated_genes) as genes:
        for rec in Bio.SeqIO.parse(genes, 'fasta'):
            gene_seqs[rec.id] = str(rec.seq).upper()
    return gene_seqs


def scan_gene_feature(features_file):
    """ Read TAB-delimited *.genes files from gene_annotations """
    features = defaultdict(dict)
    with InputStream(features_file) as stream:
        for r in select_from_tsv(stream, selected_columns=genes_feature_schema, result_structure=dict):
            if r['gene_type'] == "CDS":
                features[r['contig_id']][r['gene_id']] = r ## gnl|Prokka|
    return features


def fetch_genes_are_markers(cluster_info_fp):
    assert os.path.exists(cluster_info_fp), f"{cluster_info_fp} doesn't exit"

    list_of_markers = []
    dict_of_genes_are_markers = defaultdict(dict)

    filter_cmd = f"awk \'$8 != \"\"\'"
    with InputStream(cluster_info_fp, filter_cmd) as stream:
        for r in select_from_tsv(stream, selected_columns=["centroid_99", "marker_id"], result_structure=dict):
            dict_of_genes_are_markers[r["centroid_99"]] = r
            if r["marker_id"] not in dict_of_genes_are_markers:
                list_of_markers.append(r["marker_id"])
    return dict_of_genes_are_markers, list_of_markers


def compute_gene_boundary(features):
    """ Given list of gene ranges, generate the desired, half-open boundaries by binary search """
    gene_boundaries = defaultdict(dict)
    for contig_id in features.keys():
        features_per_contig = features[contig_id]
        # Sort gene features by starting genomic position
        dict_of_feature_tuples = {fid: (feat['start'], feat['end']) for fid, feat in features_per_contig.items()}
        feature_ranges_sorted = dict(sorted(dict_of_feature_tuples.items(), key=lambda x: x[1][0], reverse=False))
        # For non-overlapping gene ranges, linear search would report the first met range (smallest)
        # therefore, we update the gene feature files for the genes with overlapping ranges with the adjacent genes before
        prev_gene = None
        gc = 0
        gene_offsets = dict()
        for gid, grange in feature_ranges_sorted.items():
            gc += 1
            if prev_gene is None:
                prev_gene = gid
                continue
            curr_gene = gid
            curr_start = feature_ranges_sorted[curr_gene][0]
            prev_end = feature_ranges_sorted[prev_gene][1]
            if curr_start <= prev_end:
                new_curr_start = prev_end + 1
                feature_ranges_sorted[curr_gene] = (new_curr_start, feature_ranges_sorted[curr_gene][1])
                gene_offsets[curr_gene] = new_curr_start - curr_start + 1
            prev_gene = curr_gene
        # Now, we are sure the sorted feature ranges don't have overlapping anymore
        feature_ranges_flat = tuple(_ for rt in tuple(feature_ranges_sorted.values()) for _ in rt)
        # Convert ranges into half-open intervals.
        boundaries = tuple(gr + 1 if idx%2 == 1 else gr for idx, gr in enumerate(feature_ranges_flat))
        gene_boundaries[contig_id] = {"genes": list(feature_ranges_sorted.keys()), "boundaries": boundaries}
    return gene_boundaries


def cluster_short_contigs(unassigned_contigs, chunk_size):
    """ Sort the contigs in ascending order and cluster small ones. The benefit is larger contigs mostly have individual chunks"""
    sorted_contigs = {cid:clen for cid, clen in sorted(unassigned_contigs.items(), key=lambda x: x[1])}
    list_of_contigs_id = list(sorted_contigs.keys())
    list_of_contigs_length = list(sorted_contigs.values())
    istart = 0
    prev_istart = 0
    chunk_id = 0
    subset_of_contigs = defaultdict(dict)
    t = []
    while istart < len(list_of_contigs_length):
        while sum(list_of_contigs_length[prev_istart:istart+1]) <= chunk_size and istart < len(list_of_contigs_length):
            istart += 1
        chunk_id += 1
        subset_of_contigs[chunk_id] = {"chunk_id": chunk_id, "contigs_id": list_of_contigs_id[prev_istart:istart],
                                       "chunk_length": sum(list_of_contigs_length[prev_istart:istart]),
                                       "contigs_len": list_of_contigs_length[prev_istart:istart]}
        t = t + list_of_contigs_length[prev_istart:istart]
        prev_istart = istart
    assert len(t) == len(list_of_contigs_length)
    assert set(t) == set(list_of_contigs_length)
    return subset_of_contigs
