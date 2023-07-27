# Chunyu Zhao 2023-06-14
# Input: ambiguous centroid ID, gene_info.txt and genes.ffn.
# Output: clean member sequences with three intermediate files: member IDs, member raw ffn, member ambiguous ffn.

#! /usr/bin/bash

if [ $# -ne 4 ]; then
    echo "Usage: $0 CENROID.ID GENES INFO OUTPUT"
    exit 1
fi

pat="$1"
genes="$2"
info="$3"
memffn="$4"

memdir="${memffn%.ffn}"
meminfo="${memdir}.txt"

# get all the member genes ID
awk -v pat="$pat" '$2 == pat {print $1}' ${info} | awk -v pat="$pat" '$1 != pat' > $meminfo

# extract all the member genes
seqkit grep -w 0 -f ${meminfo} ${genes} > ${memffn}
