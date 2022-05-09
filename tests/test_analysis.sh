#!/bin/bash
set -e

if [ $# -ne 1 ]; then
    echo "Usage: $0 NUMCORES"
    exit 1
fi

num_cores=$1

basedir=`pwd`
testdir="${basedir}/tests"

outdir="${testdir}/midas2_output"
rm -rf ${outdir}
mkdir -p ${outdir}

midas_outdir="${outdir}/single_sample"
merge_midas_outdir="${outdir}/across_samples"

midas_dbname="uhgg"
midas_db="${outdir}/midasdb_${midas_dbname}"

logs_dir="${outdir}/logs"
mkdir -p "${logs_dir}"

samples_fp="${outdir}/samples.txt"
pool_fp="${outdir}/samples_list.tsv"

rm -rf ${samples_fp}
rm -rf ${pool_fp}

ls "${testdir}/reads" | awk -F '_' '{print $1}' > ${samples_fp}

echo "MIDAS 2.0 Unit Testing Start"

echo -e "sample_name\tmidas_outdir" > ${pool_fp}
cat ${samples_fp} | awk -v OFS='\t' -v dir=$midas_outdir '{print $1, dir}' >> ${pool_fp}


echo "Testing Single-Sample Species Module"
cat ${samples_fp} | xargs -Ixx bash -c \
    "midas2 run_species --sample_name xx -1 ${testdir}/reads/xx_R1.fastq.gz \
     --num_cores ${num_cores} --midasdb_name ${midas_dbname} --midasdb_dir ${midas_db} \
    ${midas_outdir} &> ${logs_dir}/xx_species_${num_cores}.log"


echo "Testing Across-Samples Species Module"
midas2 merge_species --samples_list ${pool_fp} \
    --median_marker_coverage 0  ${merge_midas_outdir} &> ${logs_dir}/merge_species_${num_cores}.log


echo "Testing Build Pan-Genome Bowtie2 Databases"
midas2 build_bowtie2db --midasdb_name ${midas_dbname} --midasdb_dir ${midas_db} \
    --species_list 102454 --num_cores ${num_cores} \
    --bt2_indexes_name pangenomes --bt2_indexes_dir ${merge_midas_outdir}/bt2_indexes \
    &> ${logs_dir}/build_bowtie2_pan_${num_cores}.log


echo "Testing Single-Sample SNPs Module"
cat ${samples_fp} | xargs -Ixx bash -c \
    "midas2 run_snps --sample_name xx -1 ${testdir}/reads/xx_R1.fastq.gz \
    --num_cores ${num_cores} --chunk_size 500000 \
    --midasdb_name ${midas_dbname} --midasdb_dir ${midas_db} --select_threshold=-1 \
    --advanced --ignore_ambiguous \
    --select_by median_marker_coverage,unique_fraction_covered \
    --select_threshold=5,0.5 \
    ${midas_outdir} &> ${logs_dir}/xx_snps_${num_cores}.log"


echo "Testing Across-Samples SNPs Module"
midas2 merge_snps --samples_list ${pool_fp} \
    --midasdb_name ${midas_dbname} --midasdb_dir ${midas_db} \
    --advanced --num_cores ${num_cores} --chunk_size 1000000 \
    --genome_coverage 0.7 ${merge_midas_outdir} \
    &> ${logs_dir}/merge_snps_${num_cores}.log


echo "Testing Single-Sample Genes Module With Existing Bowtie Database"
head -n 2 ${samples_fp} | xargs -Ixx bash -c \
    "midas2 run_genes --sample_name xx -1 ${testdir}/reads/xx_R1.fastq.gz --num_cores ${num_cores} \
     --midasdb_name ${midas_dbname} --midasdb_dir ${midas_db} \
     --prebuilt_bowtie2_indexes ${merge_midas_outdir}/bt2_indexes/pangenomes \
     --prebuilt_bowtie2_species ${merge_midas_outdir}/bt2_indexes/pangenomes.species \
     ${midas_outdir} &> ${logs_dir}/xx_genes_${num_cores}_w_bowtie2.log"


echo "Testing Across-Samples Genes Module"
midas2 merge_genes --samples_list ${pool_fp} --midasdb_name ${midas_dbname} --midasdb_dir ${midas_db} \
     --num_cores ${num_cores} --sample_counts 2 ${merge_midas_outdir} \
     &> ${logs_dir}/merge_genes_${num_cores}.log


echo "MIDAS 2.0 Unit Testing SUCCESS"
