
.. _midas2_wiki:

##########################################
Output Files and Directory Layout
##########################################

In this section we refer to the output directory for single-sample as ``midas_output/sample_name``, and across-samples as ``midas_output``.
Input command-line options can be found at :ref:`this page<common_cli_options>`.

Single-Sample Results Layout
============================

MIDAS 2.0 analysis usually starts with species selection which selects sufficiently abundant species in each sample (command ``run_species``).
After completing this step, users can run either of two strain-level analysis: ``run_snps`` for single-sample read pileup (SNV module) or
``run_genes`` for pan-gene profiling (CNV module).

Here is an example of the results layout of all single-sample analysis in the local filesystem.

.. code-block:: shell

  Output                                 Producer        Meaning
  -------------------------------------------------------------------------------------------
  {midas_output}/{sample_name}
    |- species
       |- species_profile.tsv            run_species     Summary of species coverage
       |- markers_profile.tsv            run_species     Per species marker coverage
    |- snps
       |- snps_summary.tsv               run_snps        Summary of reads mapping to rep-genome
       |- {species}.snps.tsv.lz4         run_snps        Per species pileups
    |- genes
       |- genes_summary.tsv              run_genes       Summary of reads mapping to pan-genome
       |- {species}.genes.tsv.lz4        run_genes       Per species pan-gene coverage

   |- temp
       |- snps
          |- repgenomes.bam              run_snps        Rep-genome alignment file
          |- {species}/snps_XX.tsv.lz4
       |- genes
          |- pangenome.bam               run_genes       Pan-genome alignment file
          |- {species}/genes_XX.tsv.lz4
    |- bt2_indexes
       |- snps/repgenomes.*              run_snps        Sample-specific rep-genome database
       |- genes/pangenomes.*             run_genes       Sample-specific pan-genome database


Across-Samples Results Layout
=============================

For a collection of samples, population SNV and pan-genome CNV can be estimated using subcommands ``merge_snps`` and ``merge_genes``.

.. code-block:: shell

  Output                                             Producer        Meaning
  ---------------------------------------------------------------------------------------------------------------
  {midas_output}
    |- species
      |- species_prevalence.tsv                      merge_species   Per species summary statistics across samples
      |- species/species_read_counts.tsv             merge_species   Species-by-sample reads counts matrix
      |- species/species_coverage.tsv                merge_species   Species-by-sample marker coverage matrix
      |- species/species_rel_abundance.tsv           merge_species   Species-by-sample relative abundance matrix
    |- snps
      |- snps_summary.tsv                            merge_snps      Alignment summary statistics per sample
      |- {species}/{species}.snps_info.tsv.lz4       merge_snps      Per species SNV metadata
      |- {species}/{species}.snps_freqs.tsv.lz4      merge_snps      Per species site-by-sample MAF matrix
      |- {species}/{species}.snps_depth.tsv.lz4      merge_snps      Per species site-by-sample reads depth matrix
    |-genes
      |- genes_summary.tsv                           merge_genes     Alignment summary statistics per sample
      |- {species}/{species}.genes_presabs.tsv.lz4   merge_genes     Per species gene-by-sample pre-abs matrix
      |- {species}/{species}.genes_copynum.tsv.lz4   merge_genes     Per species gene-by-sample copy number matrix
      |- {species}/{species}.genes_depth.tsv.lz4     merge_genes     Per species gene-by-sample reads depth matrix



.. _db_layout:

MIDAS Reference Database Layout
===============================

To meet the challenge of increased number of available genome sequences,
MIDAS 2.0 implemented a new database infrastructure, geared to run on `AWS Batch <https://aws.amazon.com/batch/>`_
and `S3 <https://aws.amazon.com/s3/>`_, to achieve `elastic scaling <https://github.com/czbiohub/pairani/wiki>`_
for building MIDAS 2.0 reference databases.

To be specific, the MIDAS 2.0 reference database construction step can be executed in AWS using hundreds
of r5d.24xlarge instances over a period of a couple of days, depositing built products in S3.
For example, it took ~$80,000 and a week to build the species pan-genome for all 47,894 species of GTDB r202.


Table of Content
----------------

The new database infrastructure reads in a table of contents (TOC) file, containing genome-to-species assignment
and a choice of representative genome for each species cluster.
One TOC file (``genomes.tsv``) per MIDAS 2.0 reference database. The TOC file has four columns,
among which ``genome_is_representative`` specify whether the ``genome`` is the representative genome
for the corresponding ``species``. Only one ``representative`` per ``species``.

.. csv-table::
  :align: left

    genome,species,representative,genome_is_representative
    GUT_GENOME138501,104351,GUT_GENOME269084,0
    GUT_GENOME269084,104351,GUT_GENOME269084,1

By default, MIDAS 2.0 inherits the representative genome assignments from published prokaryotic genome databases.
Inspired by the importance of selecting proper reference genome for accurate template-based SNP calling,
this new infrastructure empowers user the flexibility to dynamically re-assign the representative genomes,
simply by modifying the ``genomes.tsv`` file accordingly.


Microbial Genome Collections
----------------------------

Unified Human Gastrointestinal Genome (UHGG)
++++++++++++++++++++++++++++++++++++++++++++
A collection of 286,997 genomes assembled from metagenomes, isolates and single cells from human stool samples
has been clustered into 4,644 gut-only species in `UHGG 1.0 catalogues <http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/>`_.
The collection of all the UHGG genomes were mirrored in a `S3 bucket <s3://jason.shi-bucket/IGGdb2.0/clean_set/>`_,
which serves as the input to the database construction.
`Six-digit numeric species ids <s3://jason.shi-bucket/IGGdb2.0/alt_species_ids.tsv>`_ were arbitrarily assigned.
Instead of species name, these ``species_id`` are used as species identifier in all the reports generated by MIDAS 2.0.

Genome Taxonomy Database (GTDB)
+++++++++++++++++++++++++++++++++

`GTDB R06-RS202 <https://gtdb.ecogenomic.org/stats/r202>`_ contains 45,555 bacterial and 2,339 archaeal species clusters
spanning 258,406 genomes, released on April 27th, 2021. The genome members for each species cluster is
specified in the `sp_clusters_r202.tsv <https://data.ace.uq.edu.au/public/gtdb/data/releases/release202/202.0/auxillary_files/sp_clusters_r202.tsv>`_,
upon which order six-digit numeric species ids are assigned.
GTDB only provided the sequences of the representative genomes, and we downloaded all the genomes from
NCBI genomes repository using `genome_updater <https://github.com/pirovc/genome_updater>`_.


Target Layout and Construction
------------------------------

MIDAS 2.0 reference database (MIDASDB) primarily consist of three parts: rep-genome databases, pan-genome databases, and universal single copy genes (SGC) marker database.
The target layout of any MIDASDB follow the same relative structure, based on the root directory of the database.
The following toy example demonstrates the major steps to construct the MIDASDB and the target layout using
a collection of two genomes (``genome1`` and ``genome2``) from one species cluster ``species1``.

**TODO: insert image**

Inputs
++++++

The input collection of genomes need to be organized in the format as ``cleaned_genomes/<species>/<genome>/<genome>.fna``.
And the table of content ``genomes.tsv`` file needs to be generated accordingly,
with randomly assigned six-digit ``species_id``, to replace the species name.
The ``genome`` name can be kept as it is.

.. csv-table::
  :align: left

  genome,species,representative,genome_is_representative
  genome1,100001,genome2,0
  genome2,100001,genome2,1


Rep-Genome Database
+++++++++++++++++++

The genome annotation for all the genomes were done by `Prokka <https://github.com/tseemann/prokka>`_,
and the annotated genes were kept under the directory of ``genes_annotations/<species>/<genome>``.
The rep-genome databases for the SNPs module analysis only included the gene annotations and sequences for the representative genomes, as specified in the TOC.

.. code-block:: shell

  gene_annotations/100001/genome2/genome2.fna.lz4
  gene_annotations/100001/genome2/genome2.ffn.lz4
  gene_annotations/100001/genome2/genome2.genes.lz4


SCG Marker Database
+++++++++++++++++++

Marker genes are defined as universal, single-copy gene families.
MIDAS 2.0 uses a subset (15) of the `PhyEco gene families <https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0077033>`_.
The pre-computed HMM model of this set of 15 single copy genes (SCGs) are available at:

.. code-block:: shell

  s3://microbiome-pollardlab/uhgg_v1/marker_gene_models/phyeco/marker_genes.hmm.lz4
  s3://microbiome-pollardlab/uhgg_v1/marker_gene_models/phyeco/marker_genes.mapping_cutoffs.lz4

For each annotated genome, the homologs of 15 SCGs were identified with ``hmmsearch``,
as well as the mapping of gene id to corresponding marker gene id,
under the directory of ``marker_genes/phyeco/temp/<species>/<genome>``.

.. code-block:: shell

    marker_genes/phyeco/temp/100001/genome2/genome2.markers.fa
    marker_genes/phyeco/temp/100001/genome2/genome2.markers.map

For all the representative genomes, the identified marker genes were concatenated into monolithic ``marker_genes.fa``,
from which ``hs-blastn`` index would be constructed. The indexed ``marker_genes.fa`` serves as the SCG marker databases.

.. code-block:: shell

    marker_genes/phyeco/marker_genes.fa
    marker_genes/phyeco/marker_genes.fa.sa
    marker_genes/phyeco/marker_genes.fa.bwt
    marker_genes/phyeco/marker_genes.fa.sequence


Pan-Genome Database
+++++++++++++++++++

Species-level pan-genome refers to the set of non-redundant genes that represent the genetic diversity within one species cluster.

In order to construct the pan-genome database for each species, the first step if to concatenate the annotated genes
from its all genome members into ``pangenomes/100001/genes.ffn``.
The second step, which is also the most time-consuming step, is to cluster the concatenated genes based on 99% percent identity (PID)
using `vsearch <https://github.com/torognes/vsearch>`_.
Each cluster was represented by the gene at its center - centroid gene (``centroids.99.ffn``).
The ``centroid.99`` genes were further on clustered to 95, 90, ..., PID, respectively, and the mapping relationships were listed in ``centroid_info.txt``.
The top level ``centroids.ffn`` file represents the 99 percent identity clusters, and serves as the species pan-genome databases.

Reads are aligned to the pan-genome databases to determine the gene content of strains in a sample (``run_genes`` command),
and reads can optionally aggregated into gene clusters at any of the lower clustering thresholds across samples (``merge_genes`` command).

.. code-block:: shell

    pangenomes/100001/centroids.ffn
    pangenomes/100001/centroid_info.txt
