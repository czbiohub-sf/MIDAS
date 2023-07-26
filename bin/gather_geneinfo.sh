#! /usr/bin/bash

set -e

if [ $# -ne 2 ]; then
    echo "Usage: $0 FILENAME OUTPUT"
    exit 1
fi

input="$1"
output="$2"


awk -F'\t' '{
    if ($1 ~ /S/) {
        print $9, $9
    } else if ($1 ~ /H/) {
        print $9, $10
    }
}' $input > $output
