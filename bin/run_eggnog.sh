#!/bin/bash
set -e

if [ $# -ne 3 ]; then
    echo "Usage: $0 MIDASDBDIR NUMCORES SPECIES"
    exit 1
fi

midasdb_dir="$1"
num_cores="$2"
species_id="$3"

midasdb_name="newdb"

echo "Processing: $species_id"
ffn="${midasdb_dir}/pangenomes/${species_id}/centroids.ffn"

outdir="${midasdb_dir}/pangenomes_annotations/02_eggnog/${species_id}"
mkdir -p ${outdir}

base_dir="./eggnog_$species_id"
mkdir -p $base_dir/temp
mkdir -p $base_dir/scratch

EGGNOG_DATA_DIR=/pollard/data/projects/mwas-projects/dbs/eggnog_data

emapper.py \
  -i ${ffn} --itype CDS \
  -m diamond --sensmode more-sensitive \
  --data_dir ${EGGNOG_DATA_DIR} \
  --output ${species_id} --override \
  --dbmem --pfam_realign realign \
  --cpu 8 \
  --temp_dir $base_dir/temp --output_dir ${outdir} \
  --scratch_dir $base_dir/scratch

echo "${species_id} DONE"
rm -rf ${base_dir}

