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


PAN_GENE_INFO_SCHEMA = {
    "gene_id": str,
    "centroid_99": str,
    "centroid_95": str,
    "centroid_90": str,
    "centroid_85": str,
    "centroid_80": str,
    "centroid_75": str,
}


PAN_GENE_LENGTH_SCHEMA = {
    "gene_id": str,
    "genome_id": str,
    "gene_length": int,
}


CLUSTER_INFO_SCHEMA = {
    "centroid_99": str,
    "centroid_95": str,
    "centroid_90": str,
    "centroid_85": str,
    "centroid_80": str,
    "centroid_75": str,
    "centroid_99_length": int,
    "marker_id": str,
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
    "total_covered_marker": int,
    "unique_covered_marker": int,
    "ambiguous_covered_marker": int,
    "total_marker_counts": int,
    "unique_fraction_covered": float,
    "total_marker_length": int,
}


species_merge_schema = {
    "species_id": str,
    "marker_read_counts": int,
    "median_marker_coverage": float,
    "marker_coverage": float,
    "marker_relative_abundance": float,
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


species_prevalence_schema = {
    "species_id": str,
    "median_abundance": float,
    "mean_abundance": float,
    "median_coverage": float,
    "mean_coverage": float,
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
    "mean_coverage": float,
}


snps_chunk_summary_schema = {
    "species_id": str,
    "chunk_id": int,
    "contig_id": str,
    "chunk_length": int,
    "aligned_reads": int,
    "mapped_reads": int,
    "contig_total_depth": int,
    "contig_covered_bases": int
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


genes_summary_schema = {
    "species_id": str,
    "pangenome_size": int,
    "covered_genes": int,
    "fraction_covered": float,
    "mean_coverage": float,
    "aligned_reads": int,
    "mapped_reads": int,
    "marker_coverage": float,
}


genes_chunk_summary_schema = {
    "species_id": str,
    "chunk_id": str,
    "chunk_genome_size": int,
    "chunk_num_covered_genes": int,
    "chunk_coverage": float,
    "chunk_aligned_reads": float,
    "chunk_mapped_reads": int,
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
    "mean_coverage": float,
    "fraction_covered": float,
    "copy_number": float,
}


genes_are_markers_schema = {
    "centroid_99": str,
    "marker_id": str,
    "gene_coverage": float,
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
