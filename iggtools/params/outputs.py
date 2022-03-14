# The "output" or built MIDAS DB layout in S3.
# See https://github.com/czbiohub/MIDAS2.0/wiki/MIDAS-DB#target-layout-in-s3

from midas2.params import inputs

## we need to replace the igg here as well.
igg = inputs.igg

genomes = f"{igg}/genomes.tsv.lz4"

opsdir = f"{igg}/operations"
