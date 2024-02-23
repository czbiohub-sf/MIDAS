# These are the schames for the outputs of MIDAS

## merge_species
DECIMALS = ".3f"
DECIMALS3 = ".3f"
DECIMALS6 = ".6f"


MARKER_INFO_SCHEMA = {
    "species_id": str,
    "genome_id": str,
    "gene_id": str,
    "gene_length": int,
    "marker_id": str
}


BLAST_M8_SCHEMA = {
    'query': str,
    'target': str,
    'pid': float,
    'aln': int,
    'mis': float,
    'gaps': float,
    'qstart': float,
    'qend': float,
    'tstart': float,
    'tend': float,
    'evalue': float,
    'score': float,
}


GENE_LENGTH_SCHEMA = {
    "gene_id": str,
    "genome_id": str,
    "gene_length": int,
}


PANGENOME_INFO_SCHEMA = {
    "gene_id": str,
    "centroid_99": str,
    "centroid_95": str,
    "centroid_90": str,
    "centroid_85": str,
    "centroid_80": str,
    "centroid_75": str,
    "gene_length": int,
    "marker_id": str,
}


GENE_ANNOTATED_SCHEMA = {
    "gene_is_phage": int,
    "gene_is_plasmid": int,
    "gene_is_amr": int,
    "gene_is_me": int,
}


PANGENOME_ANNOTATED_SCHEMA = {**PANGENOME_INFO_SCHEMA, **GENE_ANNOTATED_SCHEMA}


PANGENOME_CLUSTER_SCHEMA = {
    "centroid_99": str,
    "centroid_95": str,
    "centroid_90": str,
    "centroid_85": str,
    "centroid_80": str,
    "centroid_75": str,
    "centroid_99_genome_counts": int,
    "centroid_99_genome_prevalence": float,
    "centroid_99_gene_counts": int,
    "centroid_99_gene_length": int,
    "centroid_99_marker_density": float,
    "centroid_99_marker_id": str,
    "phage_ratio": float,
    "plasmid_ratio": float,
    "amr_ratio": float,
    "me_ratio": float,
    "COG_category": str,
    "EC": str,
    "KEGG_Pathway": str,
    "Description": str,
    "PFAMs": str
}


def fetch_cluster_xx_info_schema(xx):
    return {
        f"centroid_{xx}": str,
        f"centroid_{xx}_genome_counts": int,
        f"centroid_{xx}_genome_prevalence": float,
        f"centroid_{xx}_gene_counts": int,
        f"centroid_{xx}_gene_length": int,
        f"centroid_{xx}_marker_density": float,
        f"centroid_{xx}_marker_id": str,
        "phage_ratio": float,
        "plasmid_ratio": float,
        "amr_ratio": float,
        "me_ratio": float,
        "COG_category": str,
        "EC": str,
        "KEGG_Pathway": str,
        "Description": str,
        "PFAMs": str
    }


def fetch_default_genome_depth(dbtype):
    if dbtype == "species":
        DEFAULT_GENOME_DEPTH = 1.0
    if dbtype == "genes":
        DEFAULT_GENOME_DEPTH = 1.0
    if dbtype == "snps":
        DEFAULT_GENOME_DEPTH = 5.0
    return DEFAULT_GENOME_DEPTH


species_profile_schema = {
    "species_id": str,
    "marker_read_counts": int,
    "median_marker_coverage": float,
    "marker_coverage": float,
    "marker_relative_abundance": float,
    "unique_fraction_covered": float,
}


species_prevalence_schema = {
    "species_id": str,
    "median_abundance": float,
    "mean_abundance": float,
    "median_coverage": float,
    "mean_depth": float,
    "sample_counts": float,
}


snps_profile_schema = {
    "species_id": str,
    "genome_length": int,
    "covered_bases": int,
    "total_depth": int,
    "aligned_reads": int,
    "mapped_reads": int,
    "fraction_covered": float,
    "mean_depth": float,
}


snps_pileup_basic_schema = {
    "ref_id": str,
    "ref_pos": int,
    "ref_allele": str,
    "depth": int,
    "count_a": int,
    "count_c": int,
    "count_g": int,
    "count_t": int,
}


snps_pileup_schema = {
    "ref_id": str,
    "ref_pos": int,
    "ref_allele": str,
    "depth": int,
    "count_a": int,
    "count_c": int,
    "count_g": int,
    "count_t": int,
    "major_allele": str,
    "minor_allele": str,
    "major_allele_freq": float,
    "minor_allele_freq": float,
    "allele_counts": int,
}


snps_info_schema = {
    "site_id": str,
    "major_allele": str,
    "minor_allele": str,
    "sample_counts": int,
    "snp_type": str,
    "rc_A": int,
    "rc_C": int,
    "rc_G": int,
    "rc_T": int,
    "sc_A": int,
    "sc_C": int,
    "sc_G": int,
    "sc_T": int,
    "locus_type": str,
    "gene_id": str,
    "site_type": str,
    "amino_acids": str
}


def fetch_genes_depth_schema(xx):
    return {
        f"cluster_{xx}_id": str,
        "gene_length": int,
        "aligned_reads": int,
        "mapped_reads": int,
        "total_depth": int,
        "mean_depth": float,
        "copy_number": float,
        "genome_prevalence": float,
        "marker_id": str,
    }


def fetch_genes_chunk_schema(xx):
    return {
        "species_id": str,
        f"c{xx}_id": str,
        f"c{xx}_length": int,
        "aligned_reads": int,
        "mapped_reads": int,
        "total_depth": int,
        "mean_depth": float,
        "copy_number": float,
        "genome_prevalence": float,
        "marker_id": str,
    }


genes_summary_schema = {
    "species_id": str,
    "pangenome_size": int,
    "covered_genes": int,
    "fraction_covered": float,
    "aligned_reads": int,
    "mapped_reads": int,
    "mean_depth": float,
    "marker_depth": float,
}


genes_info_schema = {
    "presabs": float,
    "copynum": float,
    "depth": float,
    "reads": int,
}


samples_pool_schema = {
    "sample_name": str,
    "midas_outdir": str,
}


def fetch_schema_by_dbtype(dbtype):
    if dbtype == "species":
        schema = species_profile_schema
    if dbtype == "snps":
        schema = snps_profile_schema
    if dbtype == "genes":
        schema = genes_summary_schema
    return schema


def format_data(x, decimal=DECIMALS3):
    return format(x, decimal) if isinstance(x, float) else str(x)


genes_feature_schema = {
    "gene_id": str,
    "contig_id": str,
    "start": int,
    "end": int,
    "strand": str,
    "gene_type": str,
}


md5sum_schema = {
    "db": str,
    "file_name": str,
    "species_id": str,
    "md5sum": str,
}


species_marker_profile_schema = {
    "species_id": str,
    "marker_id": str,
    "marker_length": int,
    "gene_id": str,
    "total_reads": int,
    "total_alnbps": int,
    "coverage": float,
    "uniq_reads": int,
    "ambi_reads": int,
    "uniq_alnbps": int,
    "ambi_alnbps": int,
}



COLS_GENOMAD = ['gene_id', 'contig_id', 'start', 'end', 'strand',
    'gene_type', 'contig_length', 'start_genomad', 'end_genomad',
    'gene', 'annotation_conjscan', 'annotation_amr',
    'annotation_accessions', 'annotation_description']


COLS_MEFINDER = ['gene_id', 'contig_id', 'start', 'end', 'strand',
    'gene_type', 'contig_length', 'start_mefinder', 'end_mefinder',
    'mge_no', 'prediction', 'name', 'type', 'synonyms']


COLS_RESFINDER = ['gene_id', 'contig_id', 'start', 'end', 'strand',
    'gene_type', 'contig_length', 'start_resfinder', 'end_resfinder',
    'resistance_gene', 'phenotype', 'accession_no']


COLS_EGGNOG = ['#query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs',
       'max_annot_lvl', 'COG_category', 'Description', 'Preferred_name', 'GOs',
       'EC', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 'KEGG_Reaction',
       'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy', 'BiGG_Reaction', 'PFAMs',
       'contig_id', 'start', 'end', 'strand', 'gene_type', 'contig_length']
