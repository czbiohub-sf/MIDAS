.. _snv_module:

##########################################
Module: Single Nucleotide Variant Analysis
##########################################

The Single Nucleotide Variant (SNV) module has two commands:

#. Single-sample allele tallying with the ``run_snps`` command;
#. SNV calling across all the samples, with ``merge_snps`` command.

The first step can be run in parallel. We assume users have already completed
the :ref:`species module<species_module>` and have
a species profile (e.g. ``midas2_output/sample1/species/species_profile.tsv``)
ready for each sample.
Alternatively, advanced users can pass a :ref:`pre-built representative genome index<build_custom_genome_index>`
ready for single-sample SNV analysis.


.. contents::
   :depth: 3


Single-Sample Analysis
======================

Conceptually, a typical invocation of the ``run_snps`` command proceeds by

#.  selecting sample-specific abundant species based on statistics calculated
    in the :ref:`species module <species_module>`;
#.  compiling these representative genomes and building a sample-specific
    bowtie2 index;
#.  mapping reads to this index with bowtie2;
#.  outputting a read mapping summary and pileup result for each representative
    genome (i.e. each species).

MIDAS2 purposely **holds off any filtering** or identification of variant
positions until the subsequent cross-sample analysis.

(In this document, we continue to use the :ref:`data<example_data>` from the
Quickstart as an example.)

.. code-block:: shell

  midas2 run_snps \
    --sample_name sample1 \
    -1 reads/sample1_R1.fastq.gz \
    --midasdb_name uhgg \
    --midasdb_dir my_midasdb_uhgg \
    --select_by median_marker_coverage,unique_fraction_covered \
    --select_threshold=2,0.5 \
    --num_cores 8 \
    midas2_output

See the documentation on :ref:`selecting abundant species<abundant_species_selection>`
for more information about the ``--select_by`` and ``--select_threshold`` flags.

.. tip::

   This step can be parallelized over samples (e.g. using shell background
   processes).

.. note::

  In MIDAS2, ``run_snps`` can automatically download
  the reference genomes for the selected species.

.. warning::

   (Race condition) If starting multiple calls to ``run_snps``
   simultaneously, be sure that reference genomes have already been
   :ref:`downloaded<init_db>`.
   Otherwise multiple redundant downloads may be started.


Cross-Samples Analysis
======================

After running all samples individually in this way, users can then
compute population SNVs across samples using the ``merge_snps`` command.

.. code-block:: shell

    midas2 merge_snps \
      --samples_list list_of_samples.tsv \
      --midasdb_name uhgg \
      --midasdb_dir my_midasdb_uhgg \
      --num_cores 8 \
      midas2_output/merge


Key Outputs
===========

Single-Sample
-------------

Unlike the :ref:`species <species_module>` and :ref:`CNV <cnv_module>` modules,
the single-sample outputs from the SNV module are less interesting than the
merged results (at least with the default mode).

..
    TODO: Link the merged section

Users may, however, find several files useful.

A summary of read alignment and pileups for each of the genomes included in the
(usually sample-specific) bowtie2 index is reported in
``midas2_output/sample1/snps/snps_summary.tsv``.

.. csv-table::
  :align: left

  *species_id*,*genome_length*,*covered_bases*,*total_depth*,*aligned_reads*,*mapped_reads*,*fraction_covered*,*mean_coverage*
  102506,5339468,2373275,8045342,468667,224553,0.444,3.390
  102337,2749621,2566404,47723458,1479479,1010530, 0.933,18.595

Where each columns has the following meaning:

.. code-block:: text

    species_id:       six-digit species id
    genome_length:    genome length
    covered_bases:    number of bases covered by at least one post-filtered reads
    total_depth:      total read depth across all covered_bases
    aligned_reads:    total read counts across covered_bases before post-alignment filter
    mapped_reads:     total read counts across covered_bases after post-alignment filter
    fraction_covered: fraction of covered_bases (aka horizontal genome coverage)
    mean_coverage:    mean read depth across all covered_bases (aka vertical genome coverage)


For each sample and species---e.g. here sample1 and species 102506
(*E. coli*)---the per-species read pileup is found in
``midas2_output/sample1/snps/102506.snps.tsv.lz4``.
Positions are filtered to only sites in the reference genome covered by at
least two reads.

.. note::
    Large output files are compressed with `LZ4 <http://lz4.github.io/lz4/>`_ to minimize storage requirements.

..
    TODO: Link to some LZ4 docs.

When uncompressed, the contents of this file should look like the following CSV:

.. csv-table::
  :align: left

  *ref_id*,*ref_pos*,*ref_allele*,*depth*,*count_a*,*count_c*,*count_g*,*count_t*
  gnl|Prokka|UHGG144544_1,881435,T,11,0,0,0,11
  gnl|Prokka|UHGG144544_1,881436,T,13,0,5,0,8
  gnl|Prokka|UHGG144544_1,881437,T,12,0,6,0,6

Where the columns have the following meaning:

.. code-block:: text

    ref_id:     scaffold/contig id
    ref_pos:    reference position
    ref_allele: reference nucleotide
    depth:      number of post-filtered reads
    count_a:    post-filtered read counts of A allele
    count_c:    post-filtered read counts of C allele
    count_g:    post-filtered read counts of G allele
    count_t:    post-filtered read counts of T allele

..
    TODO: Explain what the filtering is? What does post-filtered mean?


Across-Samples
--------------

A number of outputs result from the multiple samples population SNV analysis.

A merged pileup summary is found in ``midas2_output/merge/snps/snps_summary.tsv``.

.. csv-table::
    :align: left

    *sample_name*,*species_id*,*genome_length*,*covered_bases*,*total_depth*,*aligned_reads*,*mapped_reads*,*fraction_covered*,*mean_coverage*
    sample1,100122,2560878,2108551,10782066,248700,207047,0.823,5.113
    sample2,100122,2560878,2300193,39263110,1180505,820736,0.898,17.069

The reported columns from ``genome_length`` to ``mean_coverage`` are the same as from
the single-sample SNV summary.


For each species, information about SNVs identified across samples is written
to ``midas2_output/merge/snps/102506.snps_info.tsv.lz4``.

.. csv-table::
  :align: left

    *site_id*,*major_allele*,*minor_allele*,*sample_counts*,*snp_type*,*rc_A*,*rc_C*,*rc_G*,*rc_T*,*sc_A*,*sc_C*,*sc_G*,*sc_T*,*locus_type*,*gene_id*,*site_type*,*amino_acids*
    gnl|Prokka|UHGG000587_14|34360|A,A,C,2,bi,26,10,0,0,2,2,0,0,CDS,UHGG000587_02083,4D,"T\,T\,T\,T"
    gnl|Prokka|UHGG000587_11|83994|T,G,T,2,bi,0,0,11,45,0,0,2,2,IGR,None,None,None

..
    (Software) Using CSV for this output that we KNOW includes ',' characters
    in the last field seems like a mistake. Wouldn't TSV be better?

Where columns have the following meaning:

.. code-block:: text

    site_id:       unique site id, composed of ref_id|ref_pos|ref_allele
    major_allele:  most common/prevalent allele in metagenomes
    minor_allele:  second most common/prevalent allele in metagenomes
    sample_counts: number of relevant samples where metagenomes is found
    snp_type:      the number of alleles observed at site (mono,bi,tri,quad)
    rc_A:          accumulated read counts of A allele in metagenomes
    rc_C:          accumulated read counts of C allele in metagenomes
    rc_G:          accumulated read counts of G allele in metagenomes
    rc_T:          accumulated read counts of T allele in metagenomes
    sc_A:          accumulated sample counts of A allele in metagenomes
    sc_C:          accumulated sample counts of C allele in metagenomes
    sc_G:          accumulated sample counts of G allele in metagenomes
    sc_T:          accumulated sample counts of T allele in metagenomes
    locus_type:    CDS (site in coding gene), RNA (site in non-coding gene), IGR (site in intergenic region)
    gene_id:       gene identified if locus type is CDS, or RNA
    site_type:     indicates degeneracy: 1D, 2D, 3D, 4D
    amino_acids:   amino acids encoded by 4 possible alleles


A site-by-sample minor allele frequency matrix is written to
``midas2_output/merge/snps/102506.snps_freq.tsv.lz4``.

.. csv-table::
  :align: left

  *site_id*,*sample1*,*sample2*
  gnl|Prokka|UHGG000587_11|83994|T,0.692,0.837
  gnl|Prokka|UHGG000587_14|34360|A,0.300,0.269

..
    Is this statistic minor / (major + minor) or minor / total?
    Is the base in the *site_id* label the major or minor allele?
    ...Or maybe the reference genome allele?

A site-by-sample read depth matrix is written to
``midas2_output/merge/snps/102506.snps_freq.tsv.lz4``.

.. note::
    This table only accounts for the alleles matching the population major
    and/or minor allele. Other bases are dropped.

.. csv-table::
  :align: left

  *site_id*,*sample1*,*sample2*
  gnl|Prokka|UHGG000587_11|83994|T,13,43
  gnl|Prokka|UHGG000587_14|34360|A,10,26


Advanced SNV Calling
====================

Single-Sample Post-alignment Filter
-----------------------------------

Users can adjust post-alignment filters via the following command-line options (default values indicated):

- ``--mapq >= 10``: discard read alignment with alignment quality < 10
- ``--mapid >= 0.94``: discard read alignment with alignment identity < 0.94
- ``--aln_readq >= 20``: discard read alignment with mean quality < 20
- ``--aln_cov >= 0.75``: discard read alignment with alignment coverage < 0.75
- ``--aln_baseq >= 30``: discard bases with quality < 30
- ``--paired_only``: only recruit properly aligned read pairs for post-alignment filter and pileup
- ``--fragment_length 5000``: maximum fragment length for paired-end alignment. Incorrect fragment length would affect the number of proper-aligned read pairs


Single-Sample Advanced SNV Calling
----------------------------------

In recognition of the need for single-sample variant calling,
we provided ``--advanced`` option to users for single-sample variant calling for all the species in the rep-genome index
with ``run_snps`` command.

In the ``--advanced`` mode, per-species pileup results will also report major allele and minor allele
for all the genomic sites covered by at least two post-filtered reads,
upon which custom variant calling filter can be applied by the users.
Users are advised to use the setting ``--ignore_ambiguous`` to avoid falsely
calling major/minor alleles for sites with tied read counts.

.. code-block:: shell

    midas2 run_snps
      --sample_name sample1 \
      -1 reads/sample1_R1.fastq.gz \
      --midasdb_name uhgg \
      --midasdb_dir my_midasdb_uhgg \
      --select_by median_marker_coverage,unique_fraction_covered \
      --select_threshold=2,0.5 \
      --fragment_length 2000 --paired_only \
      --advanced --ignore_ambiguous \
      --num_cores 8
      midas2_output


Expected Output
^^^^^^^^^^^^^^^

In the ``--advanced`` mode, per-species pileup results will include five additional columns of the major/minor allele for all the covered genomic sites.

.. csv-table::
  :align: left

    *ref_id*,*ref_pos*,*ref_allele*,*depth*,*count_a*,*count_c*,*count_g*,*count_t*,*major_allele*,*minor_allele*,*major_allele_freq*,*minor_allele_freq*,*allele_counts*
    gnl|Prokka|UHGG144544_1,881435,T,11,0,0,0,11,T,T,1.000,0.000,1
    gnl|Prokka|UHGG144544_1,881436,T,13,0,5,0,8,T,C,0.615,0.385,2
    gnl|Prokka|UHGG144544_1,881437,T,12,0,6,0,6,C,T,0.500,0.500,2

-   ``major_allele``: the allele with the most read counts
-   ``minor_allele``: the allele with the 2nd most read counts; same with major_allele if only one allele is observed
-   ``major_allele_freq``: allele frequency of ``major_allele``
-   ``minor_allele_freq``: allele frequency of ``minor_allele``; 0.0 if only one allele is observed
-   ``allele_counts``: number of alleles observed


Adjust Population SNV Filters
-----------------------------

Advanced users can refer to :ref:`this page<population_snv_calling>` for understanding the compute of population SNV.
The species, sample, and site filters for the across-samples SNV calling can be customized with command-line options. For example,

-   We can select species with ``horizontal coverage > 40%``, ``vertical coverage > 3X`` and present in more than 2 relevant samples:

.. code-block:: shell

    --genome_coverage 0.4 --genome_depth 3 --sample_counts 2

-   We can apply the following site selections: only consider site with ``read depth >= 5``, and ``read depth <= 3 * genome_depth``, and the minimal allele frequency to call an allele present is 0.05.

.. code-block:: shell

    --site_depth 5 --site_ratio 3 --snp_maf 0.05

-   We can only report populations SNV meeting the following criteria: bi-allelic, common population SNV (present in more than 80% of the population) from the protein coding genes based on accumulated sample counts.

.. code-block:: shell

    --snp_type bi --snv_type common --site_prev 0.8 --locus_type CDS --snp_pooled_method prevalence

Now we can put all the above-mentioned filters in one `merge_snps` command:

.. code-block:: shell

    midas2 merge_snps
      --samples_list list_of_samples.tsv \
      --midasdb_name uhgg \
      --midasdb_dir my_midasdb_uhgg \
      --genome_coverage 0.4 --genome_depth 3 --sample_counts 2 \
      --site_depth 5 --site_ratio 3 --snp_maf 0.05 \
      --snp_type bi --snv_type common --site_prev 0.8 --locus_type CDS --snp_pooled_method prevalence \
      --num_cores 8 \
      midas2_output/merge
