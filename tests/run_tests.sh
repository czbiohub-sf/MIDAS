#!/bin/bash

set -e
set -x

set -o pipefail

midas_outdir="my_midas_output"
logs_dir="logs"
merge_midas_outdir="merged_midas_output"
num_cores=4

mkdir -p ${logs_dir}


echo "test midas_run_species"
cat samples.txt | xargs -Ixx bash -c "python -m iggtools midas_run_species --sample_name xx -1 reads/xx_R1.fastq.gz --num_cores ${num_cores} --debug ${midas_outdir}_default &> ${logs_dir}/xx_species_default.log"


echo "test midas_run_snps default"
cat samples.txt | xargs -Ixx bash -c "python -m iggtools midas_run_snps --sample_name xx -1 reads/xx_R1.fastq.gz --num_cores ${num_cores} --debug --marker_depth 1.0 ${midas_outdir}_default &> ${logs_dir}/xx_snps_default.log"


echo "test midas_run_genes default"
cat samples.txt | xargs -Ixx bash -c "python -m iggtools midas_run_genes --sample_name xx -1 reads/xx_R1.fastq.gz --num_cores ${num_cores} --debug --marker_depth 1.0 ${midas_outdir}_default &> ${logs_dir}/xx_genes_default.log"


echo "test midas_merge_species"
python -m iggtools midas_merge_species --samples_list samples_list.tsv ${merge_midas_outdir} &> ${logs_dir}/merge_species.log


echo "test build_bowtie2: select species by prevalence"
python -m iggtools build_bowtie2_indexes --midas_iggdb ${merge_midas_outdir} --species_profile ${merge_midas_outdir}/species/species_prevalence.tsv --select_by sample_counts --select_threshold 2 --num_cores ${num_cores} --bt2_indexes_dir ${merge_midas_outdir}/bt2_indexes --debug &> ${logs_dir}/build_bowtie2.log


echo "test midas_merge_snps default"
python -m iggtools midas_merge_snps --samples_list samples_list.tsv --num_cores ${num_cores} ${merge_midas_outdir} &>  ${logs_dir}/merge_snps.log


echo "test midas_merge_genes default"
python -m iggtools midas_merge_genes --samples_list samples_list.tsv --num_cores ${num_cores} ${merge_midas_outdir} &> ${logs_dir}/merge_genes.log


echo "test midas_run_snps with prebuilt bowtie indexes"
cat samples.txt | xargs -Ixx bash -c "python -m iggtools midas_run_species --sample_name xx -1 reads/xx_R1.fastq.gz --num_cores ${num_cores} --debug ${midas_outdir}_w_bowtie2 &> ${logs_dir}/xx_species_w_bowtie2.log"
cat samples.txt | xargs -Ixx bash -c "python -m iggtools midas_run_snps --sample_name xx -1 reads/xx_R1.fastq.gz --num_cores ${num_cores} --debug --marker_depth 1.0 --prebuilt_bowtie2_indexes ${merge_midas_outdir}/bt2_indexes/repgenomes --prebuilt_bowtie2_species ${merge_midas_outdir}/bt2_indexes/repgenomes.species ${midas_outdir}_w_bowtie2 &> ${logs_dir}/xx_snps_w_bowtie2.log"
