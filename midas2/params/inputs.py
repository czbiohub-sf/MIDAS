# These are the "raw" inputs to the MIDAS DB construction subcommands.
# See https://github.com/czbiohub/MIDAS2.0/wiki/MIDAS-DB#inputs

igg = "s3://microbiome-igg"
MIDASDB_DICT = {
    "uhgg": f"{igg}/2.0",
    "gtdb": f"{igg}/gtdb",
    "testdb": f"{igg}/testdb", #<--- dummy db
}
MIDASDB_VERSION = {
    "uhgg": f"version 1.0",
    "gtdb": f"version r202",
}
MIDASDB_NAMES = list(MIDASDB_DICT.keys())


# Check out https://github.com/czbiohub/MIDAS2.0/wiki/MIDAS-DB#marker-genes-identification
marker_set = "phyeco"
MARKER_FILE_EXTS = ["fa", "fa.bwt", "fa.header", "fa.sa", "fa.sequence", "map"]
hmmsearch_max_evalue = 1e-5
hmmsearch_min_cov = 0.00


## UHGG import genomes only.
uhggv1 = "s3://jason.shi-bucket/IGGdb2.0"
genomes2species = f"{uhggv1}/genomes2species.tab"
alt_species_ids = f"{uhggv1}/alt_species_ids.tsv"
uhgg_genomes = f"{uhggv1}/clean_set"
