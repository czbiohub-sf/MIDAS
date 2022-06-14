.. _species_module:

#########################
Module: Species Selection
#########################

Reference-based metagenotyping depends crucially on the choice of reference sequences.
Microbiome data usually contains hundreds of species in one sample,
and an ideal reference database is both representative and
comprehensive in terms of the abundant species in the sample. A badly chosen reference may suffer
both from ambiguous mapping of reads to two or more sequences or spurious cross-mapping to
incorrect sequences. Therefore, a typical MIDAS 2.0 workflow starts with a species selection step,
which filters the MIDASDB to sufficiently abundant species in each particular
sample.


.. contents::
   :depth: 3


Single-Sample Analysis
======================

MIDAS 2.0 estimates species coverage by profiling the coverage of universal,
single copy, taxonomic marker genes (SCGs, 15 per species), to quickly
determine which species are abundant in the sample.


.. warning::

  This is designed *only* to select species with sufficient coverage in each
  sample. It is not intended to quantify species abundance.


(In this document, we continue to use the :ref:`data<example_data>` from the Quickstart as an example.)

.. code-block:: shell

  midas2 run_species \
      --sample_name sample1 \
      -1 reads/sample1_R1.fastq.gz \
      --midasdb_name uhgg \
      --midasdb_dir my_midasdb_uhgg \
      --num_cores 4 \
      midas2_output

.. tip::

   Single-sample analysis step can be parallelized over samples (e.g. xargs)

.. note::

  The first time ``run_species`` is used, MIDAS 2.0 will automatically download
  the marker gene database.

.. warning::

   (Race condition) If starting multiple calls to ``run_species``
   simultaneously, be sure that the marker gene database has already been
   :ref:`downloaded<init_db>`.
   Otherwise multiple, redundant downloads may be started.


Cross-Sample Merging
=====================

We now run the ``merge_species`` command to merge the single-sample species
profiling results for the samples listed our
:ref:`samples_list<prepare_sample_list>`.

.. code-block:: shell

  midas2 merge_species \
    --samples_list list_of_samples.tsv \
    --min_cov 2 \
    midas2_output/merge


The ``--min_cov`` flag defines the minimum ``median_marker_coverage`` for
estimating species prevalence, which is output as the ``sample_counts``
statistic. See below.


Key Outputs
===========

Single-Sample
-------------

For each sample, the primary output of the ``run_species`` command is (e.g.)
``midas2_output/sample1/species/species_profile.tsv``
This file describes the
coverage of each species' marker genes in the sample.
Species are sorted in decreasing order of ``median_marker_coverage``.
Only species with more than two marker genes covered with more than two reads
(a very low bar) are reported in this file.

.. csv-table::
  :align: left

  *species_id*,*marker_read_counts*,*median_marker_coverage*,*marker_coverage*,*marker_relative_abundance*,*unique_fraction_covered*
  102337,4110,28.48,28.91,0.30,1.00
  102506,734,4.98,4.98,0.05,0.93

Where the columns have the following meaning:

.. code-block:: text

    species_id:                 six-digit species id
    marker_read_counts:         total mapped read counts
    median_marker_coverage:     median coverage of the 15 SCGs
    marker_coverage:            mean coverage of the 15 SCGs
    marker_relative_abundance:  computed based on ``marker_coverage``
    unique_fraction_covered:    the fraction of uniquely mapped SCGs genes


Downstream commands (``run_snps`` and ``run_genes``) use the
``median_marker_coverage`` and/or ``unique_fraction_covered`` to select
sufficiently abundant species. See below.


.. _across_samples_species_profile:
Across-Samples
--------------

The primary output of the merging step is the file
``midas2_output/merge/species/species_prevalence.tsv``.

.. csv-table::
  :align: left

  *species_id*,*median_abundance*,*mean_abundance*,*median_coverage*,*mean_coverage*,*sample_counts*
  102337,0.186,0.186,16.205,16.205,2
  102506,0.035,0.035,2.967,2.967,2

Where the columns have the following meaning:

.. code-block:: text

    species_id:       six-digit species id
    median_abundance: median marker_relative_abundance across samples
    mean_abundance:   mean marker_relative_abundance across samples
    median_coverage:  median median_marker_coverge across samples
    mean_coverage:    mean median_marker_coverge across samples
    sample_counts:    number of samples with median_marker_coverge >= min_cov

..
    (Software) I don't like that min_cov is a CLI flag, but not tracked anywhere
    in the output directory.
    Users who run this merge_species command and don't know to manually track
    what value they used for min_cov will have lost
    key information about how to interpret one of the columns.
    I think this is a big problem.
    I believe users should either be entirely responsible for keeping track
    of parameters AND have full control over output files, OR MIDAS can
    control complex file outputs AND MUST fully track parameters itself.
    Currently, what happens if users run merge_species with different
    values of min_cvrg? I can't tell.
    This also seems like a perfectly reasonable thing for users to do:
    run MIDAS multiple times with different parameters.


MIDAS 2.0 also writes two species-by-sample matrices in the output
directory: ``midas2_output/merge/species``.
Median marker coverage, and unique fraction covered are written to
``midas2_output/merge/species/species_marker_median_coverage.tsv`` and
``midas2_output/merge/species/species_unique_fraction_covered.tsv``, respectively


..
    (Software) Consider reformatting these outputs so that each matrix isn't a
    separate file, but rather each columns is a measure and the
    sample-by-matrix part is "stacked" into a long format.


.. _abundant_species_selection:

Species Selection in Downstream Modules
=======================================

..
    This content is shared by both SNV and CNV. We should give it its own page
    and link to it from the two modules.


In a standard SNV/CNV workflow, only sufficiently abundant species in the
restricted species profile will be included to build representative genome
(rep-genome) or pan-genome index and further to be genotyped. By default,
both the ``run_snv`` and ``run_cnv`` commands perform a species selection step.
Both commands therefore assume that ``run_species`` has already been
carried out for each sample.

Two flags, ``--select_by`` and ``--select_threshold``, determine which species are selected:

- ``--select_by`` followed by a comma separated list of column names in
  ``midas2_output/species/species_profile.tsv``
- ``--select_threshold`` followed by a comma-separated list of threshold values
  for selection.


For most analyses we recommend using the combination of
``median_marker_coverage > 2X`` and ``unique_fraction_covered > 0.5``:

.. code-block:: shell

  --select_by median_marker_coverage,unique_fraction_covered --select_threshold=2,0.5


Some users may wish to genotype low abundance species and should adjust the parameters accordingly:

.. code-block:: shell

    --select_by median_marker_coverage,unique_fraction_covered --select_threshold=0,0.5


Alternatively, users can directly pick a list of species using the ``--species_list`` option.
It is worth noting that the species in the provided species list are still subject to
the ``--select_threshold`` restriction. Users can set ``--select_threshold=-1`` to
escape species selection filters based on the species profiling:

.. code-block:: shell

    --species_list 102337,102506 --select_threshold=-1


**All** the species passing the species selection filters will be genotyped.

Having finished the species selection step, we can now go to the SNV or CNV
modules, depending on the scientific aims.
