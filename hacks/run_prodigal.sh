#!/usr/bin/env bash

set -e

if [ $# -ne 1 ]; then
    echo "Usage: $0 GENOME"
    exit 1
fi

genome="$1"

## Usage
#cat genome_prodigal.txt  | xargs -IXX -n 1 -P 96 bash -c './run_prodigal.sh XX'

echo $genome

## Run prodigal
lz4 -dc /mnt/chunyu20TB/v2.0.0/genomes/${genome}.fna.lz4 | prodigal -o /mnt/chunyu20TB/v2.0.0/prodigal/${genome}.gff -a /mnt/chunyu20TB/v2.0.0/prodigal/${genome}.faa -d /mnt/chunyu20TB/v2.0.0/prodigal/${genome}.fna &> /mnt/chunyu20TB/v2.0.0/prodigal/${genome}.log && \
 lz4 -c /mnt/chunyu20TB/v2.0.0/prodigal/${genome}.gff | aws s3 cp - s3://microbiome-igg/2.0/prodigal/${genome}.gff.lz4 && \
 lz4 -c /mnt/chunyu20TB/v2.0.0/prodigal/${genome}.faa | aws s3 cp - s3://microbiome-igg/2.0/prodigal/${genome}.faa.lz4 && \
 lz4 -c /mnt/chunyu20TB/v2.0.0/prodigal/${genome}.fna | aws s3 cp - s3://microbiome-igg/2.0/prodigal/${genome}.fna.lz4 && \
 lz4 -c /mnt/chunyu20TB/v2.0.0/prodigal/${genome}.log | aws s3 cp - s3://microbiome-igg/2.0/prodigal/${genome}.log.lz4 && \
 touch /mnt/chunyu20TB/v2.0.0/flag/prodigal_${genome}_success && \
 fetchMG.pl -m extraction -v -o /mnt/chunyu20TB/v2.0.0/fetchMG/${genome}  /mnt/chunyu20TB/v2.0.0/prodigal/${genome}.faa -t 1 && \
 cd /mnt/chunyu20TB/v2.0.0/fetchMG && \
 tar -cf - ${genome} | lz4 -c - | aws s3 cp - s3://microbiome-igg/2.0/fetchMG/${genome}.tar.lz4 && \
 touch /mnt/chunyu20TB/v2.0.0/flag/fetchMG_${genome}_success

