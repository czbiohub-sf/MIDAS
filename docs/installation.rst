.. _installation:

Installation
============

MIDAS 2.0 and all its dependencies can be installed in a few ways.

From Source
++++++++++++

`Install Conda
<https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_ if
you have not already.

Users who want the most up-to-date version of the MIDAS code
can install from source (dependencies installed with Conda).

.. code-block:: shell

  $ git clone https://github.com/czbiohub/MIDAS2.0.git
  $ cd MIDAS2.0

  $ conda env create -n midas2.0 -f midas2.yml
  $ conda activate midas2.0
  $ cpanm Bio::SearchIO::hmmer --force # Temporary fix for Prokka

  $ pip install .

.. tip::

    Using the pip ``--editable`` flag here (``pip install --editable .``)
    is useful for those wishing to modify the MIDAS source code directly.


Conda
+++++++++++++

Alternatively, users can install MIDAS 2.0 and its dependencies with conda package (``midas2``):

.. code-block:: shell

  conda config --set channel_priority flexible
  conda install -c zhaoc1 -c anaconda -c bioconda -c conda-forge -c defaults midas2



Docker
++++++++++++

We also provide a pre-built Docker container.

.. code-block:: shell

  docker pull zhaoc1/MIDAS 2.0:latest
  docker run --volume "/home/ubuntu/.aws":"/root/.aws":ro --rm -it MIDAS 2.0:latest



We've also included integration tests, which can be run using the provided
script ::

  $ bash tests/test_analysis.sh 8

This will run an example analysis with 8 cores,
and will verify that all the dependencies are correctly installed
and that all analysis modules of MIDAS 2.0 can run properly.
We recommend running this after installing MIDAS 2.0 from source.
