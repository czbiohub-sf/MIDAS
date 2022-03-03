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

outdir="${testdir}/midas_output"
mkdir -p ${outdir}
midas_outdir="${outdir}/midas_run_${num_cores}"
merge_midas_outdir="${outdir}/midas_merge_${num_cores}"
midas_db="${outdir}/my_midasdb"


logs_dir="${outdir}/logs"
mkdir -p "${logs_dir}"

samples_fp="${outdir}/samples.txt"
pool_fp="${outdir}/samples_list.tsv"

rm -rf ${samples_fp}
rm -rf ${pool_fp}

ls "${testdir}/reads" | awk -F '_' '{print $1}' > ${samples_fp}


echo -e "sample_name\tmidas_outdir" > ${pool_fp}
cat ${samples_fp} | awk -v OFS='\t' -v dir=$midas_outdir '{print $1, dir}' >> ${pool_fp}


echo "test midas_run_species"
cat ${samples_fp} | xargs -Ixx bash -c "python3 -m iggtools midas_run_species --sample_name xx -1 ${testdir}/reads/xx_R1.fastq.gz --num_cores ${num_cores} --midasdb_name uhgg --midasdb_dir ${midas_db} ${midas_outdir} &> ${logs_dir}/xx_species_${num_cores}.log"

echo "test midas_merge_species"
python3 -m iggtools midas_merge_species --samples_list ${pool_fp} --marker_depth 0.5 ${merge_midas_outdir} &> ${logs_dir}/merge_species_${num_cores}.log


echo "test build_bowtie2: select species by prevalence"
python3 -m iggtools build_bowtie2db --midasdb_name uhgg --midasdb_dir ${midas_db} --species_profile ${merge_midas_outdir}/species/species_prevalence.tsv --select_by sample_counts --select_threshold 1 --num_cores ${num_cores} --bt2_indexes_dir ${merge_midas_outdir}/bt2_indexes &> ${logs_dir}/build_bowtie2_rep_${num_cores}.log
python3 -m iggtools build_bowtie2db --midasdb_name uhgg --midasdb_dir ${midas_db} --bt2_indexes_name pangenomes --species_profile ${merge_midas_outdir}/species/species_prevalence.tsv --select_by sample_counts --select_threshold 1 --num_cores ${num_cores} --bt2_indexes_dir ${merge_midas_outdir}/bt2_indexes &> ${logs_dir}/build_bowtie2_pan_${num_cores}.log


echo "test midas_run_snps with prebuilt bowtie indexes (no need to run species flow)"
cat ${samples_fp} | xargs -Ixx bash -c "python3 -m iggtools midas_run_snps --sample_name xx -1 ${testdir}/reads/xx_R1.fastq.gz --num_cores ${num_cores} --chunk_size 500000 --midasdb_name uhgg --midasdb_dir ${midas_db} --select_threshold=-1 \
  --advanced --prebuilt_bowtie2_indexes ${merge_midas_outdir}/bt2_indexes/repgenomes --prebuilt_bowtie2_species ${merge_midas_outdir}/bt2_indexes/repgenomes.species ${midas_outdir} &> ${logs_dir}/xx_snps_${num_cores}_w_bowtie2.log"


echo "test midas_merge_snps"
python3 -m iggtools midas_merge_snps --samples_list ${pool_fp} --midasdb_name uhgg --midasdb_dir ${midas_db} --advanced --num_cores ${num_cores} --chunk_size 1000000 --genome_coverage 0.6 ${merge_midas_outdir} &>  ${logs_dir}/merge_snps_${num_cores}.log


echo "test midas_run_genes"
cat ${samples_fp} | xargs -Ixx bash -c "python3 -m iggtools midas_run_genes --sample_name xx -1 ${testdir}/reads/xx_R1.fastq.gz --num_cores ${num_cores} --midasdb_name uhgg --midasdb_dir ${midas_db} --select_by unique_fraction_covered --select_threshold 0.4 ${midas_outdir}  &> ${logs_dir}/xx_genes_${num_cores}.log"


echo "test midas_merge_genes default"
python3 -m iggtools midas_merge_genes --samples_list ${pool_fp} --midasdb_name uhgg --midasdb_dir ${midas_db} --num_cores ${num_cores} --sample_counts 2 ${merge_midas_outdir} &> ${logs_dir}/merge_genes_${num_cores}.log


echo "DONE Testing MIDAS"
