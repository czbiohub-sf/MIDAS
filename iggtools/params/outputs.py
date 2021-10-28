# The "output" or built DB layout in S3.
# See https://github.com/czbiohub/iggtools/wiki#target-layout-in-s3

from iggtools.params import inputs

## we need to replace the igg here as well.
igg = inputs.igg

genomes = f"{igg}/genomes.tsv.lz4"

opsdir = f"{igg}/operations"
