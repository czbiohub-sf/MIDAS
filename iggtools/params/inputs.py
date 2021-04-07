# These are the "raw" inputs to the IGGTools database construction subcommands.
# See https://github.com/czbiohub/iggtools/wiki#inputs

iggdb2 = "s3://jason.shi-bucket/IGGdb2.0"
genomes2species = f"{iggdb2}/genomes2species.tab"
alt_species_ids = f"{iggdb2}/alt_species_ids.tsv"
uhgg_genomes = "s3://jason.shi-bucket/IGGdb2.0/clean_set"

# Check out https://github.com/czbiohub/iggtools/wiki#marker-genes
igg = "s3://microbiome-igg/2.0"
marker_set = "phyeco"
