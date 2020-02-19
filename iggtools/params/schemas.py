# These are the schames for the outputs of MIDAS

## midas_merge_species
DEFAULT_SAMPLE_DEPTH = 1.0
DECIMALS = ".3f"

species_profile_schema = {
    "species_id": str,
    "count_reads": int,
    "coverage": float,
    "relative_abundance": float
}


## midas_merge_snps
DEFAULT_SAMPLE_COUNTS = 2
DEFAULT_GENOME_DEPTH = 5.0
DEFAULT_GENOME_COVERAGE = 0.4

DEFAULT_SITE_DEPTH = 1
DEFAULT_SITE_RATIO = 2.0
DEFAULT_SITE_PREV = 0.80
DEFAULT_SITE_MAF = 0.01
DEFAULT_ALLELE_FREQ = 0.01

DEFAULT_SNP_MAF = 0
DEFAULT_SNP_TYPE = "bi"
DEFAULT_ALLELE_TYPE = "sample_counts"

DEBUG_MAX_LINES = 1000 * 1000

DECIMALS = ".3f"


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
}

## merge_midas_genes
DEFAULT_SAMPLE_DEPTH = 1.0
DEFAULT_SAMPLE_COUNTS = 1
DEFAULT_CLUSTER_ID = '95'
DEFAULT_MIN_COPY = 0.35

DECIMALS = ".3f"

genes_profile_schema = {
    "species_id": str,
    "pangenome_size": int,
    "covered_genes": int,
    "fraction_covered": float,
    "mean_coverage": float,
    "marker_coverage": float,
    "aligned_reads": int,
    "mapped_reads": int
}


genes_info_schema = {
    "presabs": float,
    "copynum": float,
    "depth": float,
    "reads": int,
}


genes_schema = {
    "gene_id": str,
    "count_reads": int,
    "coverage": float,
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
        schema = genes_profile_schema
    return schema
