# Metagenomic Intra-Species Diversity Analysis 2

[![DOI](https://zenodo.org/badge/195910808.svg)](https://zenodo.org/badge/latestdoi/195910808)

Metagenomic Intra-Species Diversity Analysis ([MIDAS](https://genome.cshlp.org/content/26/11/1612)) is an integrated pipeline for profiling strain-level genomic variations in shotgun metagenomic data. The standard MIDAS workflow harnesses a reference database of 5,926 species extracted from 30,000 genomes (MIDAS DB v1.2). MIDAS2 used the same analysis workflow as the original [MIDAS tool](https://github.com/snayfach/MIDAS), and is engineered to work with more comprehensive MIDAS Reference Databases (MIDASDBs), and to run on  collections of thousands of samples in a fast and scalable manner.

For MIDAS2, we have already built two MIDASDBs from large, public, microbial genome databases: [UHGG 1.0](https://www.nature.com/articles/s41587-020-0603-3) and [GTDB r202](https://gtdb.ecogenomic.org/). 


Publication is available in [Bioinformatics](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btac713/6793850). User manual is available at [ReadTheDocs](https://midas2.readthedocs.io/en/latest).

The performance of reads mapping based metagenotyping pipeline depends on (1) how closely related the DB reference genomes are to the strains in the samples being genotyped, and (2) post-alignment filter options, and etc. Pitfalls of genotyping microbial communities with rapidly growing genome collections can be found [here](https://www.biorxiv.org/content/10.1101/2022.06.30.498336v1).
