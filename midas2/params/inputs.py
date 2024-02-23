# These are the "raw" inputs to the MIDAS DB construction subcommands.
# See https://github.com/czbiohub/MIDAS2.0/wiki/MIDAS-DB#inputs

igg = "s3://microbiome-pollardlab"
server="https://midasdb.pollard.gladstone.org"

MIDASDB_DICT = {
    "uhgg": f"{server}/uhgg",
    "gtdb": f"{server}/gtdb",
    "newdb": f"{server}/newdb", # reserved for building new MIDASDB locally
    "s3db": f"{igg}/testdb", # reserved for building new MIDASDB and upload to S3
}
MIDASDB_NAMES = list(MIDASDB_DICT.keys())


MIDASDB_VERSION = {
    "uhgg": f"version 2.0",
    "gtdb": f"version r202",
}

MIDASDB_STATS = {
    "uhgg": {"species":4744, "genomes":289232, "version":f"version 2.0"},
    "gtdb": {"species":47893, "genomes":258405, "version":f"version r202"}
}


MD5SUM_JSON = {
    "uhgg": "cf71fe7a82f01b81cafb4944cd3a0017",
    "gtdb": "cea856e55d0dfb82bc3d20b3b2267d5c"
}

marker_set = "phyeco"
MARKER_FILE_EXTS = ["fa", "fa.bwt", "fa.header", "fa.sa", "fa.sequence", "map"]
hmmsearch_max_evalue = 1e-5
hmmsearch_min_cov = 0.00
