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
RESFIND_DATA_DIR="/pollard/data/projects/mwas-projects/dbs/resfinder_data"

echo "Processing: $species_id:$genome_id"
ffn="${midasdb_dir}/gene_annotations/${species_id}/${genome_id}/${genome_id}.fna"

outdir="${midasdb_dir}/pangenomes_annotations/01_mge/${species_id}/${genome_id}"
mkdir -p ${outdir}

scratchdir="./scratch_$genome_id"
mkdir -p $scratchdir/temp
mkdir -p $scratchdir/output

rm -rf ${outdir}/mefinder_output
mkdir -p ${outdir}/mefinder_output
mefinder find --contig ${ffn} -t ${num_cores} --temp-dir $scratchdir/temp $scratchdir/output/mefinder
cp $scratchdir/output/mefinder* ${outdir}/mefinder_output

rm -rf ${outdir}/resfinder_output
python -m resfinder -ifa ${ffn} -o ${outdir}/resfinder_output \
  -s Other -l 0.6 -t 0.8 --acquired -db_res ${RESFIND_DATA_DIR}/resfinder_db \
  -d -db_disinf ${RESFIND_DATA_DIR}/disinfinder_db
rm -rf ${scratchdir}

