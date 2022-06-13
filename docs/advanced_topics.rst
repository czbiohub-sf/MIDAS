Advanced Topics
===============

.. contents::
   :depth: 3


.. _population_snv_calling:

Population SNV Calling
**********************

In this section, we will the compute of population SNV, the chunkified implementation of pileup, as well as filters of species, samples and genomic sites in the **SNV** module.

Important Concepts
------------------

#.  **<species, relevant samples> selection**

    Population SNV analysis restricts attention to "sufficiently well" covered species in "sufficiently many" samples.

    To be specific, a given <species, sample> pair will only be kept (by default) if it has more than 40% horizontal genome coverage (``genome_coverage``)
    and 5X vertical genome coverage (``genome_depth``).
    Furthermore, only "sufficiently prevalent" species with "sufficiently many" (``sample_counts``) would be included for the population SNV analysis.
    Therefore, different species may have different lists of relevant samples.

#.  **<site, relevant samples> selection**

    For each genomic site, a sample is considered to be "relevant" if the corresponding site depth falls between the range defined by the input
    arguments ``site_depth`` and ```site_ratio * mean_genome_coverage``; otherwise it is ignored for the across-samples SNV compute.

    Therefore, different genomic sites from the same species may have different panels of "relevant samples".
    And genomic site prevalence can be computed as the ratio of the number of relevant samples for the given site over the total number of relevant samples for the given species.

#.  **relevant site**

    For each species, a site is considered to be "relevant" if the site prevalence meets the range defined by the input arguments ``snv_type`` and ```site_prev``.
    By default, common SNV with more than 90% prevalence are reported.


.. _population_snv_computation:

Population SNV Computation
--------------------------

There are three main steps to compute and report population SNV in MIDAS 2.0.

First, for each relevant genomic site, MIDAS 2.0 determines the set of alleles present across all relevant samples.
Specifically, for each allele (A, C, G, T), ``merge_snps`` command

#. tallys the sample counts (``sc_``) of relevant samples containing corresponding allele (``scA:scT``)
#. sums up the read counts (``rc`_`) of the corresponding allele across all the relevant samples (``rc_G:rc_T``).

.. csv-table::
  :align: left

    site_id,rc_A,rc_C,rc_G,rc_T,sc_A,sc_C,sc_G,sc_T
    gnl|Prokka|UHGG000587_14|34360|A,26,10,0,0,1,2,0,0

Second, population major and minor alleles for a single site can be computed based on the
accumulated read counts or sample counts across all relevant samples.
The population major allele refers to the most abundant/prevalent allele, and
the population minor allele refers to the second most prevalent/abundant allele.

For example, the population major allele of the site ``gnl|Prokka|UHGG000587_14|34360|A`` in the above example is ``A`` defined
by accumulated read counts and ``C`` defined by accumulated sample counts.

Third, MIDAS 2.0 collects and reports the sample-by-site matrix of the corresponding (1) site depth and (2)
allele frequency of the above calculated population minor allele for all the relevant samples.
In these two matrices, MIDAS 2.0 encode ``site_depth = 0`` and ``allele_frequency = -1`` with the special meaning of missing <site, sample> pair.


Chunkified Pileup Implementation
--------------------------------

Both single-sample and across-samples pileup are parallelized on the unit of chunk of sites, which is indexed by <species_id, chunk_id>.
Only when all chunks from the same species finished processing, chunk-level pileup results will merged into species-level pileup result.

This implementation makes population SNV analysis across thousands of samples possible.
To compute the population SNV for one chunk, all the pileup results of corresponding sites across all the samples need to be read into memory.
With the uses of multiple CPUs, multiple chunks can be processed at the same time.
Therefore, for large collections of samples, we recommend higher CPU counts and smaller chunk size to
optimally balance memory and I/O usage, especially for highly prevalent species.
Users can adjust the number of sites per chunk via ``chunk_size`` (default value = 1000000).
MIDAS 2.0 also has a ``robust_chunk`` option, where assigning different chunk sizes to different species based on the species prevalence.



.. _build_your_own_database:

Build Your Own MIDASDB
**********************

MIDAS 2.0 users can locally build a new MIDASDB for a custom collection of genomes.
The target layout of MIDASDB can be found at :ref:`MIDASDB Layout<db_layout>`.
This section is focused specifically on the database construction commands.


Table of Content
----------------

To start with, users need to organize the genomes in a specific format and produce the TOC ``genomes.tsv`` as described in the MIDASDB Layout.

We have prepared a toy collections of genomes at the ``tests/genomes``, and we will build the new MIDASDB under the directory of ``tests/genomes``.
There are two command-line parameters that users need to pass:

- ``--debug``: keep the local file after successfully build the database
- ``--force``: re-build the database even if already locally exists

MIDAS 2.0 reserved the ``--midasdb_name newdb`` for building custome MIDASDB, and the new MIDASDB will be built at ``--midasdb_dir``.

Rep-genome
----------

First, annotate all the genomes:

.. code-block:: shell

  midas2 annotate_genome --species all
    --midasdb_name newdb --midasdb_dir my_new_midasdb \
    --debug --force

  midas2 build_midasdb --generate_gene_feature \
    --genomes all \
    --midasdb_name newdb --midasdb_dir my_new_midasdb
    --debug --force


SCG Markers
-----------

Second, infer SCGs for all the genomes and build marker database:

.. code-block:: shell

  midas2 infer_markers --genomes all
    --midasdb_name newdb --midasdb_dir my_new_midasdb \
    --debug --force

  midas2 build_midasdb --build_markerdb \
    --midasdb_name newdb --midasdb_dir my_new_midasdb \
    --debug --force


Pan-genome
----------

Third, build species pangenomes:

.. code-block:: shell

  midas2 build_pangenome --species all \
    --midasdb_name newdb --midasdb_dir my_new_midasdb \
    --debug --force

  midas2 build_midasdb --generate_cluster_info \
    --species all \
    --midasdb_name newdb --midasdb_dir my_new_midasdb \
    --debug --force


.. _build_custom_genome_index:

Build Your Own Genome Index
***************************


MIDAS 2.0 builds sample-specific rep-genome or pan-genome index for species in the restricted species profile.
However, we recognize the needs of using one comprehensive list of species across samples in the same study.
And in this section, we will go over the steps of building one genome index a list of customized species across a given panel of samples.

We presuppose users have already completed the :ref:`across-samples species profiling<species_module>`
and have ``midas2_output/merge/species/species_prevalence.tsv`` ready for the given panel of samples.

Species Selection
-----------------

Users can select species based on the prevalence from the ``species_prevalence.tsv`` file, e.g. the list of speices that is present in at least one sample,
by customizing the ``--select_by`` and ``--selectd_threshold`` to the ``build_bowtie2db`` command.

Build Genome Index
------------------

In this section, we will keep using the :ref:`example data<example_data>` from Quickstart.

.. code-block:: shell

  midas2 build_bowtie2db \
    --midasdb_name uhgg --midasdb_dir my_midasdb_uhgg \
    --select_by sample_counts \
    --select_threshold 2 \
    --bt2_indexes_name repgenomes \
    --bt2_indexes_dir one_bt2_indexes \
    --num_cores 8

And users can locate the generated rep-genome database at ``one_bt2_indexes/repgenomes``, and the list of species in the rep-genome is at ``one_bt2_indexes/repgenomes.species``.

Use Prebuilt Genome Index
-------------------------

If taking this approach, for the single-sample SNV or CNV analysis, users can pass the pre-built rep-genome to ``run_snps`` analysis (pan-genome for ``run_genes``), as following:

.. code-block:: shell

  midas2 run_snps
    --sample_name sample1 \
    -1 reads/sample1_R1.fastq.gz \
    --midasdb_name uhgg \
    --midasdb_dir my_midasdb_uhgg \
    --prebuilt_bowtie2_indexes one_bt2_indexes/repgenomes \
    --prebuilt_bowtie2_species one_bt2_indexes/repgenomes.species \
    --select_threshold=-1 \
    --num_cores 8 \
    ${midas_output}



Developer Notes
**********************

Export Your Conda Environment
-----------------------------

.. code-block:: shell

  conda update --all
  conda clean –all
  conda env export --no-builds | grep -v "^prefix:" > midas2.updated.yml
