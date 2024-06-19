# These are the "raw" inputs to the MIDAS DB construction subcommands.
# See https://github.com/czbiohub/MIDAS2.0/wiki/MIDAS-DB#inputs

igg = "s3://microbiome-pollardlab"
server="https://midasdb.pollard.gladstone.org"

MIDASDB_DICT = {
    "uhgg": f"{server}/v3/uhgg",
    "gtdb": f"{server}/v3/gtdb",
    "localdb": f"{server}/localdb", # reserved for building new MIDASDB locally
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
    "uhgg": "c5c05ea8747fe11c7f0b08bb69dc280b",
    "gtdb": "73c657e156f520770c57ad334cb044ef"
}

marker_set = "phyeco"
MARKER_FILE_EXTS = ["fa", "fa.bwt", "fa.header", "fa.sa", "fa.sequence", "map"]
hmmsearch_max_evalue = 1e-5
hmmsearch_min_cov = 0.00
