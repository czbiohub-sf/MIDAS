# The "output" or built DB layout in S3.
# See https://github.com/czbiohub/iggtools/wiki#target-layout-in-s3

igg = "s3://microbiome-igg/2.0"
genomes = f"{igg}/genomes.tsv"
cleaned_imports = f"{igg}/cleaned_imports"
pangenomes = f"{igg}/pangenomes"
annotations = f"{igg}/prokka"
opsdir = f"{igg}/operations"
