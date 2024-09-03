# Set the base image
FROM ubuntu:22.04
ENV DEBIAN_FRONTEND noninteractive

# File Author / Maintainer
MAINTAINER Chunyu Zhao chunyu.zhao@czbiohub.org


ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"


RUN apt-get update && apt-get -y -o Dpkg::Options::="--force-confdef" -o Dpkg::Options::="--force-confold" upgrade
RUN apt-get install -y apt-utils
RUN apt-get install -y pkg-config build-essential cpanminus git make perl cmake autoconf wget
RUN apt-get install -y mdadm xfsprogs htop # for aegea
RUN apt-get install -y sudo # for s3mi
RUN apt-get install -y nano

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


RUN conda install -c anaconda python=3.9.18 numpy pandas lz4 -y
RUN conda install -c conda-forge awscli biopython ncurses time -y
RUN conda install -c bioconda bowtie2 samtools pysam prokka hs-blastn=0.0.5 vsearch gffutils pybedtools cd-hit -y

RUN pip3 install bottleneck==1.3.6 --upgrade
RUN pip3 install midasv3
RUN midas --version


# For aegea
RUN pip3 install awscli --upgrade
RUN pip3 install aegea
RUN pip3 install 'git+https://github.com/chanzuckerberg/s3mi.git@v0.8.0'
RUN aegea --version
RUN s3mi


RUN conda clean --all -y


WORKDIR /
