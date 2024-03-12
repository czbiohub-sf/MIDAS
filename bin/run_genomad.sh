#!/bin/bash
set -e

if [ $# -ne 4 ]; then
    echo "Usage: $0 MIDASDBDIR NUMCORES SPECIES GENOME"
    exit 1
fi

midasdb_dir="$1"
num_cores="$2"
genome_id="$3"
species_id="$4"

midasdb_name="newdb"

GENOMAD_DATA_DIR="/pollard/data/projects/mwas-projects/dbs/genomad_db"

echo "Processing: $species_id:$genome_id"
ffn="${midasdb_dir}/gene_annotations/${species_id}/${genome_id}/${genome_id}.fna"

outdir="${midasdb_dir}/pangenomes_annotations/01_mge/${species_id}/${genome_id}"
mkdir -p ${outdir}
rm -rf ${outdir}/genomad_output
mkdir -p ${outdir}/genomad_output

genomad end-to-end --enable-score-calibration ${ffn} ${outdir}/genomad_output ${GENOMAD_DATA_DIR} --threads ${num_cores} --cleanup

