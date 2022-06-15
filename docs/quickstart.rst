Quickstart
============

Setup
*****


Install MIDAS2
-----------------

On a Linux machine, download a copy of MIDAS2 from our GitHub repository, and
install the dependencies. We do not currently support non-Linux environments.


..
    Confirmed that this is the best install method? Seems like the most complicated one...
    I guess part of the goal here is to get all of the test reads which are saved to
    github...


.. code-block:: shell

  $ git clone https://github.com/czbiohub/MIDAS2.git
  $ cd MIDAS2

  $ conda env create -n midas2 -f midas2.yml
  $ conda activate midas2
  $ cpanm Bio::SearchIO::hmmer --force # Temporary fix for Prokka

Install MIDAS2.

.. code-block:: shell

  $ pip install .


Alternative installation procedures are also :ref:`described elsewhere<installation>`.


.. _example_data:

Example Data
------------

Assuming you're in the ``MIDAS2.0/`` directory you just ``cd``-ed to,
two single-end gzipped FASTQ files are in the folder ``tests/reads``.
Navigate to the ``tests`` directory ::

  $ cd tests


.. _init_db:

Pre-download SCG Genes
**********************

..
    I think you should delete this pre-loading step, since
    MIDAS is designed to do it automatically.
    If you intend to remove this functionality soon, but
    otherwise I think it fits the quickstart mentality to use
    as much of the automated stuff as possible.


Download the universal single copy genes for MIDAS Reference Database (MIDAS DB) of ``uhgg``
to a new subfolder called ``my_midasdb_uhgg`` ::

  $ midas2 database --init --midasdb_name uhgg --midasdb_dir my_midasdb_uhgg

..
    TODO: If I'm not mistaken, this will install the MIDASDB to MIDAS2.0/tests/my_midasdb_uhgg
    Seems like a mistake, since users will run the quickstart and then have to _redo_ the
    database download when they want to run MIDAS on different project...
    TODO: Is there a similar issue with using the _installation_ detailed above?
    Will users need to uninstall and re-install for some reason?

..
    TODO: Add links to the more completely explanations of each step
    elsewhere in the wiki.


.. _demo_midas_ourdir:

Identify Abundant Species
*************************

We'll start by searching for reads that align to single-copy, taxonomic marker
genes in order to identify abundant species in each sample.

.. code-block:: shell

  for sample_name in sample1 sample2
  do
    midas2 run_species \
        --sample_name ${sample_name} \
        -1 reads/${sample_name}_R1.fastq.gz \
        --midasdb_name uhgg \
        --midasdb_dir my_midasdb_uhgg \
        --num_cores 4 \
        midas2_output
    done


Single-nucleotide Variant Analysis
**********************************

Identify SNVs in Each Sample
----------------------------
..
    Is "SNV calling" an accurate description of what MIDAS is doing here?
    Seems more like this step is just about alignment to the reference
    genome and SNV-calling only really happens in the cross-sample analysis.

We'll next run the single-sample SNV analysis for each sample.

.. code-block:: shell

  for sample_name in sample1 sample2
  do
    midas2 run_snps \
      --sample_name ${sample_name} \
      -1 reads/${sample_name}_R1.fastq.gz \
      --midasdb_name uhgg \
      --midasdb_dir my_midasdb_uhgg \
      --num_cores 4 \
      midas2_output
  done

The pileup summary for ``sample1`` is written to
``midas2_output/sample1/snps/snps_summary.tsv``.
This file summarizes the read mapping
and pileup results for each of the abundant species determined in the previous
step.
By default, species are selected based on the filter:
``median_marker_coverage > 2``. More details about abundant species selection can
be referred :ref:`here<abundant_species_selection>`.


Compute Population SNVs across multiple samples
-----------------------------------------------

.. _prepare_sample_list:


In order to compute population SNV from multiple single-sample pileup results, we first
need to construct a tab-separated **sample manifest file**: ``list_of_samples.tsv``.

This file has a column for the ``sample_name`` and another for
``midas_output``, and is required for multi-sample analyses.

.. code-block:: shell

  echo -e "sample_name\tmidas_outdir" > list_of_samples.tsv
  ls reads | awk -F '_' '{print $1}' | awk -v OFS='\t' '{print $1, "midas2_output"}' >> list_of_samples.tsv


..
    TODO: The shell command to build this file is a bit opaque, and users
    may have other ideas about how to build it. Maybe skip the shell
    script and just provide the manifest already in ``reads/``?

We can take a look at the ``list_of_samples.tsv``:

.. code-block:: shell

  cat list_of_samples.tsv
  sample_name	midas_outdir
  sample1	midas2_output
  sample2	midas2_output


Based on this output, we can run ``merge_snps`` and MIDAS2 will know to
look at ``midas2_output/sample1/snps/snps_summary.tsv`` for the ``run_snps``
output from sample1.


Now we are ready to compute the population SNVs across these two samples:

.. code-block:: shell

  midas2 merge_snps \
    --samples_list list_of_samples.tsv \
    --midasdb_name uhgg \
    --midasdb_dir my_midasdb_uhgg \
    --genome_coverage 0.7 \
    --num_cores 4 \
    midas2_output/merge


Users may be interested in the contents of the file
``midas2_output/merge/snps_summary.tsv`` written in this step.

.. code-block:: shell

  cat midas2_output/merge/snps_summary.tsv
  sample_name	species_id	genome_length	covered_bases	total_depth	aligned_reads	mapped_reads	fraction_covered	mean_coverage
  sample1	102454	2762447	2322823	15271923	145639	131992	0.841	6.575
  sample2	102454	2762447	2322823	15270765	145639	131982	0.841	6.574


Other output files and the full output directory structure can be found at
:doc:`output`.


Copy-number Variant Analysis
**********************************

Identify CNVs in Each Sample
----------------------------

Since building bowtie2 indexes for the species pangenomes takes much longer time, we
first build the bowtie2 indexes for one species (102454) to a new subfolder ``bt2_indexes/``:

.. code-block:: shell

  midas2 build_bowtie2db \
    --midasdb_name uhgg --midasdb_dir my_midasdb_uhgg \
    --species_list 102454 \
    --bt2_indexes_name pangenomes \
    --bt2_indexes_dir bt2_indexes \
    --num_cores 4

More information about building your own bowtie2 indexes for either representative genome (repgenome)
or pangenome can referred :ref:`here<build_custom_genome_index>`.


Now we can run the single-sample CNV analysis for each sample with the existing bowtie2 indexes.
The pileup summary for ``sample1`` will be generated under the directory
``midas2_output/sample1/genes/genes_summary.tsv``.


.. code-block:: shell

  for sample_name in sample1 sample2
  do
    midas2 run_genes \
      --sample_name ${sample_name} \
      -1 reads/${sample_name}_R1.fastq.gz \
      --midasdb_name uhgg \
      --midasdb_dir my_midasdb_uhgg \
      --prebuilt_bowtie2_indexes bt2_indexes/pangenomes \
      --prebuilt_bowtie2_species bt2_indexes/pangenomes.species \
      --num_cores 4 \
      midas2_output
  done


Compile CNVs across multiple samples
------------------------------------

Same with the population SNV analysis, multi-sample CNV analysis also requires a tab-separated
:ref:`sample manifest file<prepare_sample_list>`.


We can then merge the per-sample CNV results:

.. code-block:: shell

  midas2 merge_genes \
    --samples_list list_of_samples.tsv \
    --midasdb_name uhgg \
    --midasdb_dir my_midasdb_uhgg \
    --num_cores 4 \
    midas2_output/merge


Users may be interested in the contents of the file
``midas2_output/merge/genes_summary.tsv`` written in this step.


.. code-block:: shell

  cat midas2_output/merge/genes_summary.tsv
  sample_name	species_id	pangenome_size	covered_genes	fraction_covered	mean_coverage	aligned_reads	mapped_reads	marker_coverage
  sample1	102454	129140	4004	0.031	3.495	162476	28611	3.410
  sample2	102454	129140	4199	0.033	3.603	169286	34908	3.410


Other output files and the full output directory structure can be found at
:doc:`output`.
