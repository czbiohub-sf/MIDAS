# Set the base image
FROM ubuntu:18.04
ENV DEBIAN_FRONTEND noninteractive

# File Author / Maintainer
MAINTAINER Chunyu Zhao chunyu.zhao@czbiohub.org


ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"


RUN apt-get update && apt-get -y -o Dpkg::Options::="--force-confdef" -o Dpkg::Options::="--force-confold" upgrade
RUN apt-get install -y apt-utils
RUN apt-get install -y pkg-config build-essential cpanminus git make perl cmake autoconf wget


RUN rm -rf /var/lib/apt/lists/*


RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh
RUN conda --version


RUN conda config --add channels anaconda
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge


RUN conda install -c anaconda python=3.7.9 numpy lz4 -y
RUN conda install -c conda-forge awscli biopython -y
RUN conda install -c bioconda bowtie2=2.4.5 samtools pysam prokka hs-blastn=0.0.5 vsearch gffutils perl-app-cpanminus -y
RUN conda clean --all -y


RUN cpanm Bio::SearchIO::hmmer --force


RUN pip3 install 'git+https://github.com/czbiohub/MIDAS2.0.git' --upgrade
RUN midas2 --version


WORKDIR /
