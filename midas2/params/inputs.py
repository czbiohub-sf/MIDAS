# These are the "raw" inputs to the MIDAS DB construction subcommands.
# See https://github.com/czbiohub/MIDAS2.0/wiki/MIDAS-DB#inputs

igg = "s3://microbiome-pollardlab"
MIDASDB_DICT = {
    "uhgg": f"{igg}/uhgg_v1",
    "gtdb": f"{igg}/gtdb_r202",
    "testdb": f"{igg}/testdb", # <-- dummy db
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
