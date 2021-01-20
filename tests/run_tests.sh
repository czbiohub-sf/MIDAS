#!/bin/bash

set -e
set -x

set -o pipefail


if [ $# -ne 1 ]; then
    echo "Usage: $0 NUMCORES"
    exit 1
fi

num_cores=$1

midas_outdir="midas_output_${num_cores}"
merge_midas_outdir="merged_midas_output_${num_cores}"

logs_dir="logs"
mkdir -p ${logs_dir}

samples_fp="samples.txt"
pool_fp="samples_list.tsv"

rm -rf ${samples_fp}
rm -rf ${pool_fp}

ls reads | awk -F '_' '{print $1}' > ${samples_fp}


echo -e "sample_name\tmidas_outdir" > ${pool_fp}
cat ${samples_fp} | awk -v OFS='\t' -v dir=$midas_outdir '{print $1, dir}' >> ${pool_fp}


echo "test midas_run_species"
cat ${samples_fp} | xargs -Ixx bash -c "python -m iggtools midas_run_species --sample_name xx -1 reads/xx_R1.fastq.gz --num_cores ${num_cores} --debug ${midas_outdir} &> ${logs_dir}/xx_species_${num_cores}.log"


echo "test midas_run_snps"
cat ${samples_fp} | xargs -Ixx bash -c "python -m iggtools midas_run_snps --sample_name xx -1 reads/xx_R1.fastq.gz --num_cores ${num_cores} --debug --marker_depth 1.0 ${midas_outdir} &> ${logs_dir}/xx_snps_${num_cores}.log"


echo "test midas_run_genes"
cat ${samples_fp} | xargs -Ixx bash -c "python -m iggtools midas_run_genes --sample_name xx -1 reads/xx_R1.fastq.gz --num_cores ${num_cores} --debug --marker_depth 1.0 ${midas_outdir} &> ${logs_dir}/xx_genes_${num_cores}.log"


echo "test midas_merge_species"
python -m iggtools midas_merge_species --samples_list ${pool_fp} ${merge_midas_outdir} &> ${logs_dir}/merge_species_${num_cores}.log


echo "test build_bowtie2: select species by prevalence"
python -m iggtools build_bowtie2_indexes --midas_iggdb ${merge_midas_outdir} --species_profile ${merge_midas_outdir}/species/species_prevalence.tsv --select_by sample_counts --select_threshold 2 --num_cores ${num_cores} --bt2_indexes_dir ${merge_midas_outdir}/bt2_indexes --debug &> ${logs_dir}/build_bowtie2_${num_cores}.log


echo "test midas_merge_snps default"
python -m iggtools midas_merge_snps --samples_list ${pool_fp} --num_cores ${num_cores} ${merge_midas_outdir} &>  ${logs_dir}/merge_snps_${num_cores}.log


echo "test midas_merge_genes default"
python -m iggtools midas_merge_genes --samples_list ${pool_fp} --num_cores ${num_cores} ${merge_midas_outdir} &> ${logs_dir}/merge_genes_${num_cores}.log


echo "test midas_run_snps with prebuilt bowtie indexes"
cat ${samples_fp} | xargs -Ixx bash -c "python -m iggtools midas_run_species --sample_name xx -1 reads/xx_R1.fastq.gz --num_cores ${num_cores} --debug ${midas_outdir}_w_bowtie2 &> ${logs_dir}/xx_species_${num_cores}_w_bowtie2.log"
cat ${samples_fp} | xargs -Ixx bash -c "python -m iggtools midas_run_snps --sample_name xx -1 reads/xx_R1.fastq.gz --num_cores ${num_cores} --debug --marker_depth 1.0 --prebuilt_bowtie2_indexes ${merge_midas_outdir}/bt2_indexes/repgenomes --prebuilt_bowtie2_species ${merge_midas_outdir}/bt2_indexes/repgenomes.species ${midas_outdir}_w_bowtie2 &> ${logs_dir}/xx_snps_${num_cores}_w_bowtie2.log"

echo "DONE"
