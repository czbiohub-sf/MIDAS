#!/usr/bin/env bash
set -e

if [ $# -ne 4 ]; then
    echo "Usage: $0 GENES OUTDIR PID THREADS"
    exit 1
fi

genes=$1
outdir=$2
pid=$3
threads=$4

## calculate the interger part of pid
calc(){ awk "BEGIN { print "$*" }"; }

pid_int=`calc $pid*100`

vsearch -cluster_fast ${genes} \
        -id "${pid}" -threads ${threads} \
        -centroids "${outdir}/centroids.${pid_int}.ffn" \
        -uc "${outdir}/uclust.${pid_int}.txt"

