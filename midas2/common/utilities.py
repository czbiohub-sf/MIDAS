#!/usr/bin/env python3
import os
import re
from io import StringIO
from bisect import bisect
from collections import defaultdict
import gffutils
import Bio.SeqIO
import pandas as pd
import numpy as np

from midas2.common.utils import InputStream, OutputStream, retry, select_from_tsv, tsprint, command
from midas2.params.schemas import genes_feature_schema, PANGENOME_INFO_SCHEMA, PANGENOME_CLUSTER_SCHEMA, fetch_centroid_prevalence_schema


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


def scan_cluster_info(cluster_info_path, cxx = "centroid_99"):
    return scan_info_file(cluster_info_path, PANGENOME_CLUSTER_SCHEMA, cxx)


def scan_centroid_prev(fp, xx):
    CENTROID_PREV_SCHEMA = fetch_centroid_prevalence_schema(xx)
    return scan_info_file(fp, CENTROID_PREV_SCHEMA, f"centroid_{xx}")


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


def prune_short_c99s(dict_of_c95_per_c99, dict_of_gene_length, cutoff=0.4):
    # remove c99s shorter than 40% of the length of corresponding c95
    # Moved to midasdb using pandas
    dict_of_c99s = {}
    for c95, list_of_c99s in dict_of_c95_per_c99.items():
        c95_length = dict_of_gene_length[c95]
        filtered_c99 = [c99 for c99 in list_of_c99s if dict_of_gene_length[c99] >= cutoff * c95_length]
        for _ in filtered_c99:
            dict_of_c99s[_] = c95
    return dict_of_c99s



def augment_gene_info(centroid_info, gene_to_marker, dict_of_gene_length, gene_info_file):
    """ Augment gene_info.txt with two additional columns: gene_length and marker_id """
    with OutputStream(gene_info_file) as stream:
        stream.write("\t".join(PANGENOME_INFO_SCHEMA.keys()) + "\n")
        for gene_id, r in centroid_info.items():
            gene_len = dict_of_gene_length[gene_id]
            marker_id = gene_to_marker[gene_id] if gene_id in gene_to_marker else ""
            val = [gene_id] + list(r.values()) + [gene_len, marker_id]
            stream.write("\t".join(map(str, val)) + "\n")


def compute_marker_density(genes_info, centroid_xx = "centroid_95"):
    """ For each centroid_xx, we compute the marker density across its ALL ORF members"""
    # genes_info: Pandas DataFrame
    df = genes_info[[centroid_xx, 'marker_id']]
    grouped = df.groupby(centroid_xx)
    # compute the number of ORFs assigned with a SGC marker ID
    non_nan_count = grouped['marker_id'].count()
    total_count = grouped.size()
    # marker density := the number of ORFs assigned as markers  / the total number of ORFs memebers (ratio)
    ratio = non_nan_count / total_count
    # Convert the series to a DataFrame and rename the columns
    ratio = ratio.reset_index()
    ratio.columns = [centroid_xx, f'{centroid_xx}_marker_density']
    return ratio


def fetch_marker_ids(c99_info, centroid_xx = "centroid_95"):
    # We want to make sure the marker_id and is_centroid95_marker are consistent
    mids = c99_info[c99_info[f"is_{centroid_xx}_marker"]][[centroid_xx, 'marker_id']].drop_duplicates()
    mids.columns = [centroid_xx, f'{centroid_xx}_marker_id']
    return mids


def decorate_cluster_info(genes_info, cutoff=0.4):
    c99_df = genes_info[genes_info["gene_id"] == genes_info["centroid_99"]]
    for c_xx in ['centroid_99', 'centroid_95']:
        c_xx_sgc_density = compute_marker_density(genes_info, c_xx)
        c_xx_sgc_density[f'is_{c_xx}_marker'] = c_xx_sgc_density[f'{c_xx}_marker_density'] > cutoff
        c99_df = pd.merge(c99_df, c_xx_sgc_density, how='left', on=c_xx)
    c99_df.rename(columns={'gene_length': 'centroid_99_length'}, inplace=True)
    c99_df.drop('gene_id', axis=1, inplace=True)
    # cal marker_id for centroids
    mids = fetch_marker_ids(c99_df, "centroid_95")
    c99_df = pd.merge(c99_df, mids, how='left', on=c_xx)
    c99_df.fillna('', inplace=True)
    return c99_df


def decorate_genes_info(midas_db, genes_info_fp, species_id):
    # decorate genes_info with binary gene annotation
    genes_info = pd.read_csv(genes_info_fp, sep="\t")
    for mstep in ['genomad_virus', 'genomad_plasmid', 'mefinder', 'resfinder']:
        final_file = midas_db.get_target_layout(f"panannot_{mstep}", False, species_id)
        new_df = pd.read_csv(final_file, sep='\t')

        if "genomad_virus" == mstep:
            new_col_name = "gene_is_phage"
        if "genomad_plasmid" == mstep:
            new_col_name = "gene_is_plasmid"
        if "mefinder" == mstep:
            new_col_name = "gene_is_amr"
        if "resfinder" == mstep:
            new_col_name = "gene_is_me"

        new_df[new_col_name] = 1
        new_df = new_df[['gene_id', new_col_name]].groupby('gene_id').first().reset_index()
        # append gene level MGE status
        genes_info = pd.merge(genes_info, new_df, on='gene_id', how='left')
        genes_info[new_col_name] = np.where(genes_info[new_col_name].isna(), 0, 1)
    return genes_info


def read_eggnog_csv(filename):
    # Read the file line-by-line and exclude lines starting with '##'
    with open(filename, 'r') as file:
        lines = [line for line in file if not line.startswith('##')]
    # Convert the filtered lines into a DataFrame
    df = pd.read_csv(StringIO('\n'.join(lines)), sep='\t')
    return df


def annotation_ratio_x_members(midas_db, species_id, xx='80'):
    """ Compuate the annotation ratio for centroid_xx level across its all gene members' annotation """
    genes_annotated = pd.read_csv(midas_db.get_target_layout("pangenome_genes_annot", False, species_id), sep="\t")
    eggnog = read_eggnog_csv(midas_db.get_target_layout("panannot_eggnog", False, species_id))

    target_cols = ['gene_is_phage', 'gene_is_plasmid', 'gene_is_amr', 'gene_is_me']

    centroids_xx = f'centroid_{xx}'
    cxx_df = genes_annotated.groupby([centroids_xx])[target_cols].sum() / genes_annotated.groupby([centroids_xx])[target_cols].count()
    cxx_df.columns = ['phage_ratio', 'plasmid_ratio', 'amr_ratio', 'me_ratio']
    cxx_df = cxx_df.reset_index()

    cxx_df = pd.merge(cxx_df, eggnog[['#query', 'COG_category', 'EC', 'KEGG_Pathway', 'Description', 'PFAMs']], left_on = centroids_xx, right_on = '#query', how='left')
    cxx_df = cxx_df.drop(columns=['#query'])
    return cxx_df
