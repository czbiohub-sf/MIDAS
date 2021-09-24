# Set the base image
FROM ubuntu:18.04
ENV DEBIAN_FRONTEND noninteractive

# File Author / Maintainer
MAINTAINER MICROBIOME-IGG Contributors bdimitrov@chanzuckerberg.com


# Ubuntu
RUN apt-get update && apt-get -y -o Dpkg::Options::="--force-confdef" -o Dpkg::Options::="--force-confold" upgrade
RUN apt-get install -y apt-utils
RUN apt-get install -y pkg-config bsdtar alien build-essential libbz2-dev liblz4-tool lbzip2 zlib1g-dev zip unzip liblzma-dev
RUN apt-get install -y sysstat emacs-nox autoconf gcc g++ curl wget gdebi-core git make perl cmake


# Python 3.7
RUN apt-get update && apt-get install -y software-properties-common
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt-get update && apt-get install -y python3.7 python3.7-dev python3-pip cython
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.7 2


RUN apt-get install -y bowtie2 vsearch samtools
RUN pip3 install pysam biopython numpy gffutils


# Prokka
RUN apt-get install -y libdatetime-perl libxml-simple-perl libdigest-md5-perl default-jre git
RUN apt-get install -y bioperl hmmer xz-utils cpanminus
# We have to manually delete version 2.6 installed by bioperl and install 2.12.0 for prokka.
WORKDIR /tmp
RUN apt-get remove -y ncbi-blast+
RUN wget -N ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.12.0/ncbi-blast-2.12.0+-1.x86_64.rpm
RUN alien -i ncbi-blast-2.12.0+-1.x86_64.rpm
RUN rm -rf ncbi-blast-2.12.0+-1.x86_64.rpm
# The following line differs from the official instructions in that it actually works.
RUN cpanm List::Util
WORKDIR /usr/local
RUN git clone https://github.com/tseemann/prokka.git && /usr/local/prokka/bin/prokka --setupdb
RUN ln -sf /usr/local/prokka/bin/prokka /usr/local/bin/prokka
RUN prokka --version


# HS-BLASTN
WORKDIR /usr/local
RUN git clone https://github.com/chenying2016/queries.git && \
  cd queries/hs-blastn-src && make
RUN ln -s /usr/local/queries/Linux-amd64/bin/hs-blastn /usr/local/bin/hs-blastn


# For aegea
RUN apt-get install -y mdadm xfsprogs htop
RUN pip3 install awscli --upgrade
RUN pip3 install --upgrade 'git+git://github.com/chanzuckerberg/s3mi.git'
# We need sudo to exist for some s3mi commands, even though it doesn't do anything
RUN apt-get install -y sudo


# This layer re-installs tbl2asn, which is a component of blast that expires after 6 months,
# and is needed by prokka.  Force rebuild of this layer when your jobs start failing in Prokka
# with error message indicating that tbl2asn has expired.
#RUN wget ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux.tbl2asn.gz && \
#    gunzip linux.tbl2asn.gz && \
#    mv linux.tbl2asn tbl2asn && \
#    chmod +x tbl2asn && \
#    mv tbl2asn /usr/local/prokka/binaries/linux/


# Cleanup
RUN rm -rf /tmp/*

WORKDIR /
