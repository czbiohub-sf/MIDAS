# The "output" or built DB layout in S3.
# See https://github.com/czbiohub/iggtools/wiki#target-layout-in-s3

from iggtools.params import inputs

igg = inputs.igg
genomes = f"{igg}/genomes.tsv.lz4"

cleaned_imports = f"{igg}/cleaned_imports"
pangenomes = f"{igg}/pangenomes"
annotations = f"{igg}/gene_annotations"

opsdir = f"{igg}/operations"
marker_genes = f"{igg}/marker_genes/{inputs.marker_set}"
marker_centroids = f"{marker_genes}/marker_centroids_99"
