# Chunyu Zhao 2023-06-30
# Input: replace tsv assigned centroids with ffn assigned centroids

#! /usr/bin/bash
set -e

if [ $# -ne 2 ]; then
    echo "Usage: $0 PATTERN TSV"
    exit 1
fi

cid="$1"
cdhit_tsv="$2"

wrong_cid=`awk -v pat="$cid" '$2 == pat {print $3}' $cdhit_tsv`

# replace the mis-labled centroids, with the correct one
awk -v pat="$wrong_cid" -v newpat="$cid" -F'\t' 'BEGIN {OFS = FS} $3 == pat {$3 = newpat}1' $cdhit_tsv > $cdhit_tsv.temp
mv $cdhit_tsv.temp $cdhit_tsv
