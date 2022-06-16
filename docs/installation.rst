.. _installation:

Installation
============

MIDAS2 and all its dependencies can be installed in a few ways.

Conda
+++++++++++++

`Install Conda
<https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_ if
you have not already. Users can install MIDAS2 and its dependencies with conda package (``midas2``):

.. code-block:: shell

  conda config --set channel_priority flexible
  conda install -c zhaoc1 -c anaconda -c bioconda -c conda-forge -c defaults midas2


If this installation takes a long time to solve dependcies conflicts,


From Source
++++++++++++

Alternatively, users who want the most up-to-date version of the MIDAS code
can install from source (dependencies installed with Conda).

On a Linux machins, download a copy of MIDAS2 from our GitHub repository,
and install the dependencies.
We do not currently support non-Linux environments.

.. code-block:: shell

  $ git clone https://github.com/czbiohub/MIDAS2.git
  $ cd MIDAS2

  $ conda env create -n midas2.0 -f midas2.yml
  $ conda activate midas2
  $ cpanm Bio::SearchIO::hmmer --force # Temporary fix for Prokka

  $ pip install .


.. tip::

    Using the pip ``--editable`` flag here (``pip install --editable .``)
    is useful for those wishing to modify the MIDAS source code directly.



Docker
++++++++++++

We also provide a pre-built Docker container.

.. code-block:: shell

  docker pull zhaoc1/MIDAS2:latest
  docker run --volume "/home/ubuntu/.aws":"/root/.aws":ro --rm -it MIDAS2:latest



We've also included integration tests, which can be run using the provided
script ::

  $ bash tests/test_analysis.sh 8

This will run an example analysis with 8 cores,
and will verify that all the dependencies are correctly installed
and that all analysis modules of MIDAS2 can run properly.
We recommend running this after installing MIDAS2 from source.
