#!/usr/bin/env python3
import re
from io import StringIO
from bisect import bisect
from collections import defaultdict
from functools import reduce
import gffutils
import Bio.SeqIO
import pandas as pd
import numpy as np

from midas2.common.utils import InputStream, OutputStream, retry, select_from_tsv, tsprint, command
from midas2.params.schemas import genes_feature_schema, PANGENOME_CLUSTER_SCHEMA


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
                        gid_int = int(re.search(r'\d+$', gid).group())
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
    # Translation table 11: none differences from the standard code
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


def reverse_complement(seq):
    """ Reverse complement sequence """
    return ''.join([complement(base) for base in list(seq[::-1])])


def extract_sequence_by_position(genome_seq, start, end, strand):
    seq = genome_seq[start-1 : end]
    if strand == "-":
        return reverse_complement(seq)
    return seq


def get_contig_length(fasta_file):
    clen = {}
    with InputStream(fasta_file) as file:
        for rec in Bio.SeqIO.parse(file, 'fasta'):
            clen[rec.id] = len(rec.seq)
    return clen


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
    assert len(curr_seq) % 3 == 0, f"gene {curr_gene_id} must by divisible by 3 to id codons"

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
def scan_info_file(info_file, schema, key):
    """
    Scan gene_info or cluster_info into memory.
     - schema := PANGENOME_CLUSTER_SCHEMA
     - key := gene_id, centroid_99
    """
    dict_of_info = {}
    with InputStream(info_file) as stream:
        for r in select_from_tsv(stream, selected_columns=schema, result_structure=dict):
            kid = r[key]
            if kid not in dict_of_info:
                dict_of_info[kid] = r
    return dict_of_info


def scan_cluster_info(cluster_info_path, xx="99", schema=PANGENOME_CLUSTER_SCHEMA):
    cxx = f"centroid_{xx}"
    return scan_info_file(cluster_info_path, schema, cxx)


@retry
def scan_genes(annotated_genes):
    """" Lookup of seq_id to sequence for PATRIC genes """
    gene_seqs = {}
    with InputStream(annotated_genes) as genes:
        for rec in Bio.SeqIO.parse(genes, 'fasta'):
            gene_seqs[rec.id] = str(rec.seq).upper()
    return gene_seqs


@retry
def scan_gene_feature(features_file):
    """ Read TAB-delimited *.genes files from gene_annotations """
    features = defaultdict(dict)
    with InputStream(features_file) as stream:
        for r in select_from_tsv(stream, selected_columns=genes_feature_schema, result_structure=dict):
            features[r['contig_id']][r['gene_id']] = r ## gnl|Prokka|
    return features


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
        gene_offsets = {}
        for gid, _ in feature_ranges_sorted.items():
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


def has_ambiguous_bases(sequence):
    # Check if sequence contains lower-case letters, which usually indicate soft-masked bases
    ambiguous_bases = ['N', 'X', 'n', 'x']
    return any(base in ambiguous_bases for base in sequence)


def parse_gff_to_tsv(gff3_file, genes_file):
    """ Convert GFF3 features format into genes.feature """
    command(f"rm -f {genes_file}")
    command(f"rm -f {gff3_file}.db")
    db = gffutils.create_db(gff3_file, f"{gff3_file}.db")
    with OutputStream(genes_file) as stream:
        stream.write("\t".join(["gene_id", "contig_id", "start", "end", "strand", "gene_type"]) + "\n")
        for feature in db.all_features():
            if feature.source == "prokka":
                continue
            if "ID" not in feature.attributes: #CRISPR
                continue
            seqid = feature.seqid
            start = feature.start
            stop = feature.stop
            strand = feature.strand
            gene_id = feature.attributes['ID'][0]
            locus_tag = feature.attributes['locus_tag'][0]
            assert gene_id == locus_tag
            gene_type = feature.featuretype
            stream.write("\t".join([gene_id, seqid, str(start), str(stop), strand, gene_type]) + "\n")
    return True


def extract_genomeid(s):
    # Replace everything after the right most "_" with empty string.
    index = s.rfind('_')
    return s[:index] if index != -1 else s


def compute_cxx_length(genes_info, centroid_xx = "centroid_99", func='median'):
    """ For each centroid_xx, we compute either (1) median gene length or (2) max gene length across all gene members """
    if func == "median":
        cxx_length = genes_info[[centroid_xx, 'gene_length']].groupby(centroid_xx)['gene_length'].median().reset_index()
    else:
        cxx_length = genes_info[[centroid_xx, 'gene_length']].groupby(centroid_xx)['gene_length'].max().reset_index()
    cxx_length.columns = [centroid_xx, f'{centroid_xx}_gene_length']
    return cxx_length



def compute_cxx_gene_counts(genes_info, centroid_xx="centroid_99"):
    c_xx_gene_counts = genes_info[[centroid_xx, 'gene_id']].groupby(centroid_xx).size().reset_index(name='total_gene_counts')
    c_xx_gene_counts.columns = [centroid_xx, f'{centroid_xx}_gene_counts']
    return c_xx_gene_counts


def list_cxx_coordinate_to_genome(genes_info, centroid_xx="centroid_99", qry_genome = ""):
    """ This function only locate the coordinate of centroids back to the query genome. """
    df = genes_info[['gene_id', centroid_xx]].copy()
    df['genome_id'] = df['gene_id'].apply(lambda x: extract_genomeid(x))
    cxx_lifted_to_qry = df[df['genome_id'] == qry_genome]
    cxx_lifted_to_qry = cxx_lifted_to_qry[['gene_id', centroid_xx]]
    # TODO: to identify multi-copy genes, we need to blast the lifetd pan-genes against qry genome
    return cxx_lifted_to_qry


def render_full_stacked_cxx_by_genome(genes_info, c_xx='centroid_99'):
    df = genes_info[['gene_id', c_xx]].copy()
    df['genome_id'] = df['gene_id'].apply(lambda x: extract_genomeid(x))
    centroid_xx_info = df[[c_xx, 'genome_id']].drop_duplicates()
    return centroid_xx_info


def compute_cxx_prevalence(genes_info, c_xx='centroid_99'):
    """ For each centroid_xx, we compute the operational gene family / cluster prevalence across all genomes """
    # Parse genome_id from the gene_id
    df = genes_info[['gene_id', c_xx]].copy()
    df['genome_id'] = df['gene_id'].apply(lambda x: extract_genomeid(x))
    # Subset the genes_info with custom list of genomes
    total_genomes = df['genome_id'].nunique()
    # Group_by centroid_xx and count the number of genomes with genes falling into given operational gene cluster
    centroid_xx_info = df[[c_xx, 'genome_id']].drop_duplicates()
    centroid_xx_prev = centroid_xx_info.groupby(c_xx)['genome_id'].count().reset_index(name='genome_counts')
    centroid_xx_prev['prevalence'] = centroid_xx_prev['genome_counts'] / total_genomes
    centroid_xx_prev.columns = [c_xx, f"{c_xx}_genome_counts", f"{c_xx}_genome_prevalence"]
    return centroid_xx_prev


def select_marker_by_max_counts(df, centroid_xx):
    # This is more like a remedy than a solution.
    # Edge case 1: cxx has two both marker_id and empty marker_id (NaN)
    # Edge case 2: rarely a erroneous fusion ORF could be assigned to two markers, and we keep the marker_id assignment with max gene_counts
    # Edge case 3: we drop ties
    df['max_counts'] = df.groupby(centroid_xx)['occurrence'].transform('max')
    max_occur_df = df[df['occurrence'] == df['max_counts']]
    max_occur_df = max_occur_df.drop(columns=['max_counts'])
    max_occur_df = max_occur_df[~max_occur_df.duplicated(subset=centroid_xx, keep=False)]
    return max_occur_df


def impute_cxx_marker_id(genes_info, c_xx="centroid_99", cutoff=0.5):
    """ Marker assignment for each centroid_xx is voted by all gene members """
    df = genes_info[[c_xx, 'marker_id']]
    # First compute the total_gene_counts per centroid_xx
    gene_counts = df.groupby(c_xx).size().reset_index(name='total_gene_counts')
    # Second compute the centroid_xx, marker_id occurrences, exlucding NA marker_ids
    pair_counts = df.dropna(subset=['marker_id']).groupby([c_xx, 'marker_id']).size().reset_index(name='occurrence')
    max_pair_counts = select_marker_by_max_counts(pair_counts, c_xx)
    # Third compute the ratio
    c_xx_markers = pd.merge(max_pair_counts, gene_counts, on=c_xx)
    c_xx_markers['marker_density'] = c_xx_markers['occurrence'] / c_xx_markers['total_gene_counts']
    # Fouth only keep marker assignment if (max) occurrence > cutoff
    passed_markers = c_xx_markers[c_xx_markers['marker_density'] > cutoff][[c_xx, "marker_id"]]
    # Fifth we still want to keep the low marker_density
    c_xx_markers = pd.merge(c_xx_markers[[c_xx, 'marker_density']], passed_markers, on=c_xx, how='left')
    c_xx_markers.columns = [c_xx, f'{c_xx}_marker_density', f'{c_xx}_marker_id']
    return c_xx_markers


def decorate_cxx_with_gene_info(genes_info, centroid_xx = 'centroid_99'):
    cxx_summary = genes_info.groupby(centroid_xx)['gene_length'].agg(
        mean_gene_length='mean',
        median_gene_length='median',
        std_gene_length='std',
        min_gene_length='min',
        max_gene_length='max',
    ).reset_index()
    return cxx_summary


def generate_cluster_xx_info(genes_info, c_xx = 'centroid_99', cutoff=0.5):
    # Say we start with c99, which is the original cluster_info
    c_xx_df = genes_info[genes_info["gene_id"] == genes_info[c_xx]]
    c_xx_prev = compute_cxx_prevalence(genes_info, c_xx)
    c_xx_gene_counts = compute_cxx_gene_counts(genes_info, c_xx)
    c_xx_length = compute_cxx_length(genes_info, c_xx, 'max') #<---
    c_xx_markers = impute_cxx_marker_id(genes_info, c_xx, cutoff)
    #c_xx_gene_info = decorate_cxx_with_gene_info(genes_info, c_xx)
    if c_xx == 'centroid_99':
        c_xx_df = c_xx_df.drop(['gene_id', 'gene_length', 'marker_id'], axis=1)
    else:
        c_xx_df = c_xx_df[c_xx]
    list_of_dfs = [c_xx_df, c_xx_prev, c_xx_gene_counts, c_xx_length, c_xx_markers]
    result = reduce(lambda left, right: pd.merge(left, right, on=c_xx, how='left'), list_of_dfs)
    # Replace NaN with 0 for marker_density
    result[f'{c_xx}_marker_density'] = result[f'{c_xx}_marker_density'].fillna(0)
    return result


def fetch_new_col_name(mstep, by="gene"):
    if mstep == "genomad_virus":
        new_col_name = f"{by}_is_phage"
    if mstep == "genomad_plasmid":
        new_col_name = f"{by}_is_plasmid"
    if mstep == "mefinder":
        new_col_name = f"{by}_is_amr"
    if mstep == "resfinder":
        new_col_name = f"{by}_is_me"
    return new_col_name


def decorate_genes_info_with_annot(genes_info, annot_files):
    # decorate genes_info with binary gene annotation
    for mstep, mfile in annot_files.items():
        new_col_name = fetch_new_col_name(mstep)
        annot_df = pd.read_csv(mfile, sep='\t')
        annot_df[new_col_name] = 1
        annot_df = annot_df[['gene_id', new_col_name]].groupby('gene_id').first().reset_index()
        # append gene level MGE status
        genes_info = pd.merge(genes_info, annot_df, on='gene_id', how='left')
        genes_info[new_col_name] = np.where(genes_info[new_col_name].isna(), 0, 1)
    return genes_info


def scan_eggnog(filename):
    # Read the file line-by-line and exclude lines starting with '##'
    with open(filename, 'r') as file:
        lines = [line for line in file if not line.startswith('##')]
    # Convert the filtered lines into a DataFrame
    df = pd.read_csv(StringIO('\n'.join(lines)), sep='\t')
    return df


def annotation_ratio_x_members(genes_annotated, eggnog_file, xx='75'):
    """ Compuate the annotation ratio for centroid_xx level across its all gene members' annotation """
    eggnog = scan_eggnog(eggnog_file)
    eggnog.drop_duplicates(inplace=True)

    centroids_xx = f'centroid_{xx}'
    target_cols = ['gene_is_phage', 'gene_is_plasmid', 'gene_is_amr', 'gene_is_me']
    cxx_df = genes_annotated.groupby([centroids_xx])[target_cols].sum() / genes_annotated.groupby([centroids_xx])[target_cols].count()
    cxx_df.columns = ['phage_ratio', 'plasmid_ratio', 'amr_ratio', 'me_ratio']
    cxx_df = cxx_df.reset_index()

    eggnog_cols = ['#query', 'COG_category', 'EC', 'KEGG_Pathway', 'Description', 'PFAMs']
    cxx_df = pd.merge(cxx_df, eggnog[eggnog_cols], left_on = centroids_xx, right_on = '#query', how='left')
    cxx_df = cxx_df.drop(columns=['#query'])
    return cxx_df
