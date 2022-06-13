Glossary
========

-   **Module**: i.e. one of Species, SNPs, Genes (or maybe also MIDASDB building?)

-   **Analysis**: Either the single-sample alignment-based tallying of reads belonging to species, SNVs, or genes (depending on which module) OR the merging, filtering, and summarization of these results across samples.

-   **Workflow**: The overall conceptual order of events from an implicit database and pre-requisite shotgun metagenomic samples through species selection and SNPs/genes modules, to results that will go into some sort of downstream analysis.

-   **Reference database**: The upstream UHGG, GTDB, something user supplies, etc. that is then pre-processed to make a specific...

-   **MIDASDB**: The pre-processed reference data in a format to be consumed by the MIDAS 2.0 modules, including marker genes, reference genomes, pangenomes, etc.

-   **genome Index**: Rather than bowtie2 database or some other ambiguous term

-   **Species-level pangenome**: refers to the set of non-redundant genes (centroids) clustered from all the genomes within one species cluster.
