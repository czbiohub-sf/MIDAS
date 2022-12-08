.. MIDAS2.0 documentation master file, created by
   sphinx-quickstart on Wed May 25 18:30:13 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MIDAS2
====================================


Metagenomic Intra-Species Diversity Analysis System 2 (MIDAS2) is an
integrated pipeline for profiling single nucleotide variants (SNVs) and gene
copy number variants (CNVs) in shotgun metagenomic reads. MIDAS2 implements
the same analyses as the original
`MIDAS <https://github.com/snayfach/MIDAS>`_,
but re-engineered to addresses the computational challenges presented by
increasingly large reference genome databases.

MIDAS2 was developed by `Chunyu Zhao <chunyu.zhao@czbiohub.org>`_
and Boris Dimitrov in the `Pollard Lab <https://docpollard.org/>`_ at
Chan Zuckerberg Biohub.
MIDAS2 expands on the original MIDAS developed by
`Stephen Nayfach <snayfach@gmail.com>`_.

For MIDAS2, we have already built two MIDASDBs from large, public, microbial genome databases:
`UHGG 1.0 <https://www.nature.com/articles/s41587-020-0603-3>`_ and
`GTDB r202 <https://gtdb.ecogenomic.org/>`_.


Source code is `available on GitHub
<https://github.com/czbiohub/MIDAS2.0.git>`_.

Publication is available in `Bioinformatics <https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btac713/6793850>`_.

Contents
--------

.. toctree::
    :maxdepth: 1

    quickstart
    installation
    overview
    species_module
    snv_module
    cnv_module
    download_midasdb
    output
    advanced_topics
    example_questions
    glossary
    acknowledgements
