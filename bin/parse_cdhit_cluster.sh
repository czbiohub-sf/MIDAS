#! /usr/bin/bash

set -e

if [ $# -ne 2 ]; then
    echo "Usage: $0 FILENAME OUTPUT"
    exit 1
fi

input="$1"
output="$2"


awk '/>Cluster/{
    cluster = substr($0, 10)
    next
}
{
    sub(/\.\.\. /, "\t", $0)
    sub(/>/, "", $0)
    sub(/at/, "", $0)
    sub(/,/, "", $0)
    print cluster"\t"$0
}' ${input} | awk '{$1=$1; print}' OFS="\t" > ${output}
