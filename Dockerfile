# Set the base image
FROM ubuntu:18.04
ENV DEBIAN_FRONTEND noninteractive

# File Author / Maintainer
MAINTAINER MICROBIOME-IGG Contributors bdimitrov@chanzuckerberg.com

# Ubuntu
RUN apt-get update && apt-get -y -o Dpkg::Options::="--force-confdef" -o Dpkg::Options::="--force-confold" upgrade
RUN apt-get install -y apt-utils
RUN apt-get install -y pkg-config bsdtar alien build-essential libbz2-dev liblz4-tool lbzip2 zlib1g-dev zip unzip
RUN apt-get install -y sysstat emacs-nox autoconf gcc g++ curl wget gdebi-core git make perl cmake

# python 3.7
RUN apt-get install -y software-properties-common
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt-get update && apt-get install -y python3.7
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.7 2
RUN apt-get install -y python3-pip

# aws
RUN pip3 install awscli --upgrade
RUN pip3 install --upgrade 'git+git://github.com/chanzuckerberg/s3mi.git'

# bioinformatics
RUN apt-get install -y samtools bowtie2 vsearch
RUN pip3 install pysam biopython

# Prokka
RUN apt-get install -y libdatetime-perl libxml-simple-perl libdigest-md5-perl default-jre git
RUN apt-get install -y bioperl
# We have to manually delete version 2.6 installed by bioperl and install 2.9 for prokka.
WORKDIR /tmp
RUN apt-get remove -y ncbi-blast+
RUN wget -N ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-1.x86_64.rpm
RUN alien -i ncbi-blast-2.9.0+-1.x86_64.rpm
RUN rm -rf ncbi-blast-2.9.0+-1.x86_64.rpm
RUN bash -c "yes | cpan Bio::Perl"
WORKDIR /usr/local
RUN git clone https://github.com/tseemann/prokka.git && /usr/local/prokka/bin/prokka --setupdb
RUN ln -sf /usr/local/prokka/bin/prokka /usr/local/bin/prokka
RUN prokka --version

# AWS instance setup
RUN apt-get install -y mdadm xfsprogs htop

# We need sudo to exist for some s3mi commands, even though it doesn't do anything
RUN apt-get install -y sudo

# Cleanup
RUN rm -rf /tmp/*

WORKDIR /
