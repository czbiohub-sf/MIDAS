#!/usr/bin/env bash

set -e

if [ $# -ne 1 ]; then
    echo "Usage: $0 GENOME"
    exit 1
fi

#cat genome_prodigal.txt  | xargs -IXX -n 1 -P 64 bash -c './run_hmm.sh XX'

# 20191016: already have Prodigal local, so no need to grab from S3 for now.

genome=$1

echo $genome

hmmsearch --noali --cpu 1 --domtblout marker_genes/temp/${genome}.hmmsearch /mnt/chunyu20TB/v2.0.0/phyeco.hmm prodigal/${genome}.faa
touch marker_genes/temp/${genome}.inprogress

aws s3 cp marker_genes/temp/${genome}.hmmsearch s3://microbiome-igg/2.0/marker_genes/temp/${genome}.hmmsearch
aws s3 cp marker_genes/temp/${genome}.inprogress s3://microbiome-igg/2.0/marker_genes/temp/${genome}.inprogress


# we need MIDAS in PATH and PYTHONPATH
python build_marker_db.py \
  -genome ${genome} -gene prodigal/${genome}.fna \
  -hmm marker_genes/temp/${genome}.hmmsearch \
  -hmm-cutoff /mnt/chunyu20TB/v2.0.0/phyeco.mapping_cutoffs \
  -mapfile /mnt/chunyu20TB/v2.0.0/mapfile \
  -output-fa marker_genes/temp/${genome}.phyeco.fa \
  -output-map marker_genes/temp/${genome}.phyeco.map


aws s3 cp marker_genes/temp/${genome}.phyeco.fa s3://microbiome-igg/2.0/marker_genes/temp/${genome}.phyeco.fa
aws s3 cp marker_genes/temp/${genome}.phyeco.map s3://microbiome-igg/2.0/marker_genes/temp/${genome}.phyeco.map

touch marker_genes/temp/${genome}.phyeco.success
aws s3 cp marker_genes/temp/${genome}.phyeco.success s3://microbiome-igg/2.0/marker_genes/temp/${genome}.phyeco.success
aws s3 rm s3://microbiome-igg/2.0/marker_genes/temp/${genome}.inprogress
