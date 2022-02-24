# Microbiome - Integrated Gut Genome Tools

For up to date information, please see https://github.com/czbiohub/iggtools/wiki

## Installation

```
git clone https://github.com/czbiohub/iggtools iggtools
cd iggtools
python setup.py build 
python setup.py install
bash tests/run_midas.sh 8
```

## Conda environment

Add proper channle order in the `.condarc` file, when downloading Conda for the first time.
```
conda config --add channels anaconda
conda config --add channels bioconda
conda config --add channels conda-forge
```

Create Conda environment.

```
conda env update --name=iggtools --quiet --file iggtools.yml
```

