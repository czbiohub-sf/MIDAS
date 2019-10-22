#!/usr/bin/env bash

set -e

if [ $# -ne 2 ]; then
    echo "Usage: $0 REPGENOME OUTDIR"
    exit 1
fi

# 20191018
# USAGES: grep -v "alt" alt_species_ids.tsv  | head -2 | cut -f2 | xargs -Irep bash -c "echo rep && ./pipeline_pangenome.sh rep /mnt/chunyu20TB/v2.0.0/pangenomes"

#repgenome="GUT_GENOME000001"
#output_dir="/mnt/chunyu20TB/v2.0.0/pangenome_test"

repgenome=$1
output_dir=$2

flag="${output_dir}/success/${repgenome}"
if [[ -f $flag ]] ; then
    echo "$flag exists"
    exit 0
fi

echo "NOW START for $repgenome"
## concatenate genes from all the genomes for each representative genome
mkdir -p ${output_dir}/${repgenome}

## For the given representation genomes, write gene.seq and gene.len to separate files first
mkdir -p ${output_dir}/${repgenome}/tmp_cleaned

genomes=`awk -v pat=${repgenome} '$2 == pat {print}' /mnt/chunyu20TB/v2.0.0/mapfile  | cut -f1`

echo $genomes | tr ' ' '\n' | xargs -Igg -P 2 \
 bash -c "/mnt/chunyu20TB/v2.0.0/write_genes.py -gene-file /mnt/chunyu20TB/v2.0.0/prodigal/gg.fna \
          -output-dir ${output_dir}/${repgenome}/tmp_cleaned"


echo $genomes | tr ' ' '\n' | xargs -Igg bash -c "cat ${output_dir}/${repgenome}/tmp_cleaned/gg.genes.fna" \
 > $output_dir/${repgenome}/genes.ffn
echo $genomes | tr ' ' '\n' | xargs -Igg bash -c "cat ${output_dir}/${repgenome}/tmp_cleaned/gg.genes.lens" \
 > $output_dir/${repgenome}/gene_lens.txt


## uclust_genes_99
mkdir -p ${output_dir}/${repgenome}/temp
/mnt/chunyu20TB/v2.0.0/run_uclust.sh $output_dir/${repgenome}/genes.ffn ${output_dir}/${repgenome}/temp 0.99 2


## uclust_genes
echo -e "0.95\n0.90\n0.85\n0.80\n0.75" | \
  xargs -Ixx -P 5 bash -c "/mnt/chunyu20TB/v2.0.0/run_uclust.sh \
        $output_dir/$repgenome/temp/centroids.99.ffn $output_dir/$repgenome/temp xx 1"

## store (uclusted) genes info
/mnt/chunyu20TB/v2.0.0/store_gene_info.py \
  -uclust-dir ${output_dir}/${repgenome}/temp \
  -output-file ${output_dir}/${repgenome}/gene_info.txt


## copy_centroids
cp ${output_dir}/${repgenome}/temp/centroids.99.ffn ${output_dir}/${repgenome}/centroids.ffn

## generate status
touch ${output_dir}/success/${repgenome}

rm -r ${output_dir}/${repgenome}/tmp_cleaned
aws s3 cp --recursive ${output_dir}/${repgenome} s3://microbiome-igg/2.0/pangenomes/${repgenome}
aws s3 cp ${output_dir}/success/${repgenome} s3://microbiome-igg/2.0/pangenomes/success/${repgenomes}

echo "NOW FINISH for $repgenome"
