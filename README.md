# Metagenomic Intra-Species Diversity Analysis

## MIDAS 2

[![DOI](https://zenodo.org/badge/195910808.svg)](https://zenodo.org/badge/latestdoi/195910808)

Metagenomic Intra-Species Diversity Analysis ([MIDAS](https://genome.cshlp.org/content/26/11/1612)) is an integrated pipeline for profiling strain-level genomic variations in shotgun metagenomic data. The standard MIDAS workflow harnesses a reference database of 5,926 species extracted from 30,000 genomes ([MIDAS DB v1.2](midas_db_v1.2.tar.gz)). MIDAS2 used the same analysis workflow as the original [MIDAS tool](https://github.com/snayfach/MIDAS), and is engineered to work with more comprehensive MIDAS Reference Databases (MIDASDBs), and to run on  collections of thousands of samples in a fast and scalable manner.

For MIDAS2, we have already built two MIDASDBs from large, public, microbial genome databases: [UHGG 1.0](https://www.nature.com/articles/s41587-020-0603-3) and [GTDB r202](https://gtdb.ecogenomic.org/). 


Publication is available in [Bioinformatics](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btac713/6793850). User manual is available at [ReadTheDocs](https://midas2.readthedocs.io/en/latest).

The performance of reads mapping based metagenotyping pipeline depends on (1) how closely related the DB reference genomes are to the strains in the samples being genotyped, and (2) post-alignment filter options, and etc. Pitfalls of genotyping microbial communities with rapidly growing genome collections can be found [here](https://www.biorxiv.org/content/10.1101/2022.06.30.498336v1).

Quick Installation:

```
conda create -n midas2 -c zhaoc1 -c conda-forge -c bioconda -c anaconda -c defaults midas
```

## MIDAS version 3

[![DOI](https://zenodo.org/badge/195910808.svg)](https://zenodo.org/doi/10.5281/zenodo.6637089)

MIDAS version 3, previously known as MIDAS2, features major updates to its pangenome database. These updates include a refinded curation process and a comprehensive functional annotation pipeline. MIDASDB can construct species-level pangenome databases from external reference genome collections, e.g. UHGG or GTDB, by clustering predicted genes into operational gene families (OGFs) at various average nucleotide identity (ANI) thresholds, with representative gene sequences of each OGF assigned as the centroids by vsearch.

1. MIDAS v3 made significant changes to the curation pipeline aiming to minimize the impact of fragmented gene sequences, spurious gene calls, chimeric assemblies, and redundant OGFs resulting from errors from cross-species contamination and highly fragmented MAGs.
2. Functional annotation includes a voting mechanism to assess the ratio of genes in each OGF related to phages, plasmids, mobile elements, and antimicrobial resistance, which is an improvement over common methods that relied on single centroid genes.
3. For pangenome profiling, MIDAS v3 compiles representative gene sequences from the 99% ANI level OGFs into a Bowtie2 index for alignment and quantification. It also  prunes potentially spurious singletons at 75% level or/and short OGFs. Vertical gene family coverage is calculated as the number of aligned reads over the gene length. 

&nbsp;&nbsp;&nbsp;&nbsp; The first step is to generate the pruned centroids sequences for species of interests. 

```
midas prune_centroids --midasdb_name localdb --midasdb_dir /path/to/midasdb-uhgg-v2 -t 1 --remove_singleton --species 100001 --force --debug
```

&nbsp;&nbsp;&nbsp;&nbsp;The second step is to pass the arguments to `run_genes`

```
midas run_genes --midasdb_name localdb --midasdb_dir /path/to/midasdb-uhgg-v2 --num_cores 8 --select_threshold=-1 --species_list 100723,104323,100041 --prune_centroids --remove_singleton midas_output
```


Details of these updates can be found at the provided [link](https://www.biorxiv.org/content/10.1101/2024.04.10.588779v1).

Quick Installation:

```
conda config --set channel_priority flexible
conda create -n midasv3 -c zhaoc1 -c conda-forge -c bioconda -c anaconda -c defaults midasv3=1.0.0
bash tests/test_analysis.sh 8
```
