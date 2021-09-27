#!/bin/bash
set -e

if [ $# -ne 1 ]; then
    echo "Usage: $0 NUMCORES"
    exit 1
fi

num_cores=$1

basedir=`pwd`
testdir="${basedir}/tests"
echo ${testdir}

midas_outdir="${testdir}/midas_output_${num_cores}"
merge_midas_outdir="${testdir}/merged_output_${num_cores}"
midas_db="${testdir}/midas_db"


logs_dir="${testdir}/logs"
mkdir -p "${logs_dir}"

samples_fp="${testdir}/samples.txt"
pool_fp="${testdir}/samples_list.tsv"

rm -rf ${samples_fp}
rm -rf ${pool_fp}

ls "${testdir}/reads" | awk -F '_' '{print $1}' > ${samples_fp}

echo -e "sample_name\tmidas_outdir" > ${pool_fp}
cat ${samples_fp} | awk -v OFS='\t' -v dir=$midas_outdir '{print $1, dir}' >> ${pool_fp}


echo "test midas_run_species"
cat ${samples_fp} | xargs -Ixx bash -c "python3 -m iggtools midas_run_species --sample_name xx -1 ${testdir}/reads/xx_R1.fastq.gz --num_cores ${num_cores} --midas_db ${midas_db} ${midas_outdir} &> ${logs_dir}/xx_species_${num_cores}.log"


echo "test midas_run_snps"
cat ${samples_fp} | xargs -Ixx bash -c "python3 -m iggtools midas_run_snps --sample_name xx -1 ${testdir}/reads/xx_R1.fastq.gz --num_cores ${num_cores} --midas_db ${midas_db} --select_by median_marker_coverage --select_threshold 0.5 ${midas_outdir} &> ${logs_dir}/xx_snps_${num_cores}.log"


echo "test midas_run_genes"
cat ${samples_fp} | xargs -Ixx bash -c "python3 -m iggtools midas_run_genes --sample_name xx -1 ${testdir}/reads/xx_R1.fastq.gz --num_cores ${num_cores} --midas_db ${midas_db} --select_by median_marker_coverage --select_threshold 0.5 --cache ${midas_outdir}  &> ${logs_dir}/xx_genes_${num_cores}.log"


echo "test midas_merge_species"
python3 -m iggtools midas_merge_species --samples_list ${pool_fp} --marker_depth 0.5 ${merge_midas_outdir} &> ${logs_dir}/merge_species_${num_cores}.log


echo "test build_bowtie2: select species by prevalence"
python3 -m iggtools build_bowtie2_indexes --midas_db ${midas_db} --species_profile ${merge_midas_outdir}/species/species_prevalence.tsv --select_by sample_counts --select_threshold 2 --num_cores ${num_cores} --bt2_indexes_dir ${merge_midas_outdir}/bt2_indexes &> ${logs_dir}/build_bowtie2_${num_cores}.log


echo "test midas_merge_snps default"
python3 -m iggtools midas_merge_snps --samples_list ${pool_fp} --midas_db ${midas_db} --num_cores ${num_cores} --genome_coverage 0.6 ${merge_midas_outdir} &>  ${logs_dir}/merge_snps_${num_cores}.log


echo "test midas_merge_genes default"
python3 -m iggtools midas_merge_genes --samples_list ${pool_fp} --midas_db ${midas_db} --num_cores ${num_cores} --sample_counts 2 ${merge_midas_outdir} &> ${logs_dir}/merge_genes_${num_cores}.log


echo "test midas_run_snps with prebuilt bowtie indexes"
cat ${samples_fp} | xargs -Ixx bash -c "python3 -m iggtools midas_run_species --sample_name xx -1 ${testdir}/reads/xx_R1.fastq.gz --num_cores ${num_cores} --midas_db ${midas_db} ${midas_outdir}_w_bowtie2 &> ${logs_dir}/xx_species_${num_cores}_w_bowtie2.log"
cat ${samples_fp} | xargs -Ixx bash -c "python3 -m iggtools midas_run_snps --sample_name xx -1 ${testdir}/reads/xx_R1.fastq.gz --num_cores ${num_cores} --midas_db ${midas_db} --select_by median_marker_coverage --select_threshold 0.5 --prebuilt_bowtie2_indexes ${merge_midas_outdir}/bt2_indexes/repgenomes --prebuilt_bowtie2_species ${merge_midas_outdir}/bt2_indexes/repgenomes.species ${midas_outdir}_w_bowtie2 &> ${logs_dir}/xx_snps_${num_cores}_w_bowtie2.log"


echo "DONE"
