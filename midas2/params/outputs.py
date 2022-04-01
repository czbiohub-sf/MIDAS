# The "output" or built MIDAS DB layout in S3.
# See https://github.com/czbiohub/MIDAS2.0/wiki/MIDAS-DB#target-layout-in-s3

from midas2.params.inputs import MIDASDB_DICT

def genomes(midasdb_name="uhgg"):
    igg = MIDASDB_DICT[midasdb_name]
    return f"{igg}/genomes.tsv.lz4"

def get_opsdir(midasdb_name="uhgg"):
    igg = MIDASDB_DICT[midasdb_name]
    return f"{igg}/operations"

opsdir = get_opsdir("uhgg")
