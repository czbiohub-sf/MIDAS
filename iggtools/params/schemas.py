# These are the schames for the outputs of MIDAS

## midas_merge_species
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


PAN_GENE_INFO_SCHEMA = {
    "gene_id": str,
    "centroid_99": str,
    "centroid_95": str,
    "centroid_90": str,
    "centroid_85": str,
    "centroid_80": str,
    "centroid_75": str,
}


def fetch_default_genome_depth(dbtype):
    if dbtype == "species":
        DEFAULT_GENOME_DEPTH = 1.0
    if dbtype == "genes":
        DEFAULT_GENOME_DEPTH = 1.0
        DEFAULT_SAMPLE_DEPTH = 1.0
    if dbtype == "snps":
        DEFAULT_GENOME_DEPTH = 5.0
    return DEFAULT_GENOME_DEPTH


species_profile_schema = {
    "species_id": str,
    "read_counts": int,
    "coverage": float,
    "relative_abundance": float
}


species_prevalence_schema = {
    "species_id": str,
    "median_abundance": float,
    "mean_abundance": float,
    "median_coverage": float,
    "mean_coverage": float,
    "sample_counts": float,
}

## midas_merge_snps
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


snps_profile_schema = {
    "species_id": str,
    "genome_length": int,
    "covered_bases": int,
    "total_depth": int,
    "aligned_reads": int,
    "mapped_reads": int,
    "fraction_covered": float,
    "mean_coverage": float,
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


## merge_midas_genes
DEFAULT_SAMPLE_DEPTH = 1.0
DEFAULT_SAMPLE_COUNTS = 1
DEFAULT_CLUSTER_ID = '95'
DEFAULT_MIN_COPY = 0.35

DECIMALS = ".3f"

genes_summary_schema = {
    "species_id": str,
    "pangenome_size": int,
    "covered_genes": int,
    "fraction_covered": float,
    "mean_coverage": float,
    "aligned_reads": int,
    "mapped_reads": int,
    "marker_depth": float
}


genes_info_schema = {
    "presabs": float,
    "copynum": float,
    "depth": float,
    "reads": int,
}


genes_coverage_schema = {
    "gene_id": str,
    "gene_length": int,
    "aligned_reads": int,
    "mapped_reads": int,
    "total_depth": float,
    "copy_number": float,
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


def format_data(x, decimal = DECIMALS3):
    return format(x, decimal) if isinstance(x, float) else str(x)


genes_feature_schema = {
    "gene_id": str,
    "contig_id": str,
    "start": int,
    "end": int,
    "strand": str,
    "gene_type": str,
}
