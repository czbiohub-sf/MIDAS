
.. _cnv_module:

####################################
Module: Copy Number Variant Analysis
####################################

Similar to the SNV module, the Copy Number Variant (CNV) module has two commands:

#. Single-sample, pan-genome copy number quantification with the ``run_genes`` command;
#. Merging these results across all samples with the ``merge_genes`` command.

The first step can be run in parallel.
We assume users have already completed the :ref:`species module<species_module>`
and have a species profile (e.g. ``midas2_output/sample1/species/species_profile.tsv``)
ready for each sample.
Alternatively, advanced users can pass a pre-built pangenome bowtie2 index.


.. contents::
   :depth: 3


Single-Sample Analysis
======================

Conceptually, a typical invocation of the ``run_genes`` command proceeds by

#.  selecting sample-specific abundant species based on statistics calculated
    in the :ref:`species module <species_module>`;
#.  compiling species pan-genomes and building a sample-specific
    bowtie2 index;
#.  mapping reads to this index with bowtie2;
#.  outputting a read mapping summary and copy number profile for each
    pan-genome (i.e. each species).

For each species, copy numbers for each gene are estimated as follows:

#.  For each gene in the species pangenome, read alignment metrics,
    e.g ``mapped_reads`` and ``mean_coverage``, are computed from the bowtie2
    results;
#.  species-level coverage is calculated as the median read coverage of the 15
    universal single-copy taxonomic marker genes.
#.  gene copy number (and presence/absence) is estimated by normalizing the
    gene coverage by the species coverage.

..
    Why "**e.g.** mapped_reads`` ...". Are there others? What are they? How are they used?

(In this document, we continue to use the :ref:`data<example_data>` from the
Quickstart as an example.)

.. code-block:: shell

  midas2 run_genes \
    --sample_name sample1 \
    -1 reads/sample1_R1.fastq.gz \
    --midasdb_name uhgg \
    --midasdb_dir my_midasdb_uhgg \
    --select_by median_marker_coverage,unique_fraction_covered \
    --select_threshold=2,0.5 \
    --num_cores 8 \
    midas2_output

See the documentation on `selecting abundant species <abundant_species_selection>`_
for more information about the ``--select_by`` and ``--select_threshold`` flags.

.. tip::

   This step can be parallelized over samples (e.g. using shell background
   processes).

.. note::

  In MIDAS2 ``run_genes`` can automatically download
  gene collections for the selected species.

.. warning::

   (Race condition) If starting multiple calls to ``run_genes``
   simultaneously, be sure that the gene collections have already been
   :ref:`downloaded<init_db>`.
   Otherwise multiple redundant downloads may be started.
   TODO: Link to the preload instructions here.


Cross-Samples Merging
=====================

Having run the single-sample CNV analysis for all the samples listed in the
``list_of_samples.tsv``, users next can merge the results and product a summary
using the ``merge_genes`` command.

.. code-block:: shell

    midas2 merge_genes \
      --samples_list list_of_samples.tsv \
      --midasdb_name uhgg \
      --midasdb_dir my_midasdb_uhgg \
      --num_cores 8 \
      midas2_output/merge


Key Outputs
===========

Single-Sample
-------------

For each sample (e.g. here sample1)
a summary of read alignment and CNV calling across all analyzed species
is written to ``midas2_output/sample1/genes/genes_summary.tsv``.

.. csv-table::
  :align: left

   *species_id*,*pangenome_size*,*covered_genes*,*fraction_covered*,*mean_coverage*,*aligned_reads*,*mapped_reads*,*marker_coverage*
   102337,15578,4468,0.287,16.213,1650361,450353,20.213
   102506,731186,4733, 0.006,3.803,681335,37272,2.140


Where each columns has the following meaning:

.. code-block:: text

    species_id:       six-digit species id
    pangenome_size:   number of centroids (non-redundant genes) in the species pangenome
    covered_genes:    number of centroids covered with at least one post-filtered read
    fraction_covered: fraction of covered_genes over pangenome_size
    mean_coverage:    average read depth across covered_genes
    aligned_reads:    total number of aligned reads before post-alignment filter
    mapped_reads:     total number of aligned reads after post-alignment filter
    marker_coverage:  average read depth across 15 universal SCGs in the species pangenome


Copy-number estimates are written to
``midas2_output/sample1/genes/102506.genes.tsv.lz4``
and include all genes covered by at least two reads.

.. note::
    Large output files are compressed with `LZ4 <http://lz4.github.io/lz4/>`_ to minimize storage requirements.

.. csv-table::
  :align: left

   *gene_id*,*gene_length*,*aligned_reads*,*mapped_reads*,*mean_coverage*,*fraction_covered*,*copy_number*
   UHGG143901_00483,555,14,6,2.961538,0.234234,1.384035
   UHGG143901_03589,384,103,57,32.840708,0.294271,15.347667
   UHGG143902_04031,207,9,2,1.737500,0.386473,0.811997

Where columns have the following meaning:

.. code-block:: text

    gene_id:          centroid id in the species pan-genome
    gene_length:      gene length
    aligned_reads:    number of aligned reads to gene_id before post-alignment filter
    mapped_reads:     number of aligned reads to gene_id after post-alignment filter
    mean_coverage:    average read depth of gene_id based on mapped_reads (total_gene_depth / covered_bases)
    fraction_covered: proportion of the gene_id covered by at least one read (covered_bases / gene_length)
    copy_number:      estimated copy number of gene_id based on mapped_reads (mean_coverage / median_marker_coverage)

Across-Samples
--------------

Merging across samples produces several outputs.

CNV results merged across samples are written to
``midas2_output/merge/genes/genes_summary.tsv``

.. csv-table::
  :align: left

  *sample_name*,*species_id*,*pangenome_size*,*covered_genes*,*fraction_covered*,*mean_coverage*,*aligned_reads*,*mapped_reads*,*marker_coverage*
  sample1,100122,29165,2535,0.087,4.723,263395,53006,1.435
  sample2,100122,9165,3212,0.110,16.095,1447684,263878,10.713

Besides ``sample_name``, which indexes the entries, the other
columns (``pangenome_size`` through ``marker_coverage``) are the same as in the
per-sample genes summary output.

For each species, a matrix of gene-by-sample copy-number
estimates---here species 102506 (*E. coli*)---are written to
``midas2_output/merge/genes/102506.genes_copynum.tsv.lz4``.

.. csv-table::
  :align: left

  *gene_id*,*sample1*,*sample2*
  UHGG000587_00401,33.969154,19.891455
  UHGG000587_01162,5.703398,2.821237
  UHGG000587_00962,2.370930,0.289325

Similarly, a presence absence matrix is written to
``midas2_output/merge/genes/102506.genes_preabs.tsv.lz4``.

.. csv-table::
  :align: left

   *gene_id*,*sample1*,*sample2*
   UHGG000587_00401,1,1
   UHGG000587_01162,1,1
   UHGG000587_00962,1,0


Raw vertical coverage data is reported in the same matrix form in
``midas2_output/merge/genes/102506.genes_depth.tsv.lz4``.


.. csv-table::
  :align: left

  *gene_id*,*sample1*,*sample2*
  UHGG000587_00401,48.747945,213.090622
  UHGG000587_01162,8.184746,30.222978
  UHGG000587_00962,3.402439,3.099448


Advanced CNV Calling
====================

Single-Sample Post-alignment Filter
-----------------------------------

Users can adjust post-alignment quality filter parameters via the command-line options (default vlaues indicated):

-  ``--mapq >= 2``: reads aligned to more than one genomic locations equally well are discarded (MAPQ=0,1)
-  ``--mapid >= 0.94``: discard read alignment with alignment identity < 0.94
-  ``--aln_readq >= 20``: discard read alignment with mean quality < 20
-  ``--aln_cov >= 0.75``: discard read alignment with alignment coverage < 0.75


Adjust Population CNV Filters
-----------------------------

The default ``merge_genes`` results are reported for pan-genes clustered at 95% identity (``cluster_pid``).
It further quantify the presence/absence for pan-genes by comparing the ``copy_number`` with the
user-defined minimal gene copy number (``min_copy``).
``cluster_pid`` and ``min_copy`` can be customized with the following command-line options:

- ``--genome_depth``: filter out species with ``mean_coverage`` < 1X.
- ``--min_copy``: genes with ``copy_number`` >= 0.35 are classified as present.
- ``--cluster_pid``: gene CNV results can be reported at various clustering cutoffs {75, 80, 85, 90, 95, 99}.
