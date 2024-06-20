#!/bin/bash
set -e
set -x

if [ $# -ne 1 ]; then
    echo "Usage: $0 NUMCORES"
    exit 1
fi

num_cores=$1

basedir=`pwd`
testdir="${basedir}/tests"

outdir="${testdir}/midas_output"
rm -rf ${outdir}
mkdir -p ${outdir}

midas_outdir="${outdir}/single_sample"
merge_midas_outdir="${outdir}/across_samples"

midas_dbname="gtdb"
midas_db="${outdir}/midasdb_${midas_dbname}"

logs_dir="${outdir}/logs"
mkdir -p "${logs_dir}"

samples_fp="${outdir}/samples.txt"
pool_fp="${outdir}/samples_list.tsv"

rm -rf ${samples_fp}
rm -rf ${pool_fp}

ls "${testdir}/reads" | awk -F '_' '{print $1}' > ${samples_fp}

echo "MIDASv3 Unit Testing Start"

echo -e "sample_name\tmidas_outdir" > ${pool_fp}
cat ${samples_fp} | awk -v OFS='\t' -v dir=$midas_outdir '{print $1, dir}' >> ${pool_fp}


echo "Testing Single-Sample Species Module"
cat ${samples_fp} | xargs -Ixx bash -c \
    "midas run_species --sample_name xx -1 ${testdir}/reads/xx_R1.fastq.gz \
     --num_cores ${num_cores} --midasdb_name ${midas_dbname} --midasdb_dir ${midas_db} \
    ${midas_outdir} &> ${logs_dir}/xx_species_${num_cores}.log"


echo "Testing Across-Samples Species Module"
midas merge_species --samples_list ${pool_fp} --min_cov 2  ${merge_midas_outdir} &> ${logs_dir}/merge_species_${num_cores}.log


echo "Testing Single-Sample SNV Module"
cat ${samples_fp} | xargs -Ixx bash -c \
    "midas run_snps --sample_name xx -1 ${testdir}/reads/xx_R1.fastq.gz \
    --num_cores ${num_cores} --chunk_size 1000000 \
    --midasdb_name ${midas_dbname} --midasdb_dir ${midas_db} \
    --ignore_ambiguous \
    --select_by median_marker_coverage,unique_fraction_covered \
    --select_threshold=5,0.5 \
    ${midas_outdir} &> ${logs_dir}/xx_snps_${num_cores}.log"


echo "Testing Across-Samples SNV Module"
midas merge_snps --samples_list ${pool_fp} \
    --midasdb_name ${midas_dbname} --midasdb_dir ${midas_db} \
    --num_cores ${num_cores} --chunk_size 100000 \
    --genome_coverage 0.7 ${merge_midas_outdir} \
    &> ${logs_dir}/merge_snps_${num_cores}.log


echo "Testing Build Pan-Genome Bowtie2 Databases"
midas build_bowtie2db --midasdb_name ${midas_dbname} --midasdb_dir ${midas_db} \
    --species_profile  ${merge_midas_outdir}/species/species_prevalence.tsv \
    --select_by sample_counts --select_threshold 1 --num_cores ${num_cores} \
    --bt2_indexes_name pangenomes --bt2_indexes_dir ${merge_midas_outdir}/bt2_indexes \
    &> ${logs_dir}/build_bowtie2_pan_${num_cores}.log


echo "Testing Single-Sample CNV Module With Existing Bowtie Database"
head -n 2 ${samples_fp} | xargs -Ixx bash -c \
    "midas run_genes --sample_name xx -1 ${testdir}/reads/xx_R1.fastq.gz --num_cores ${num_cores} \
     --midasdb_name ${midas_dbname} --midasdb_dir ${midas_db} --select_threshold=-1 \
     --prebuilt_bowtie2_indexes ${merge_midas_outdir}/bt2_indexes/pangenomes \
     --prebuilt_bowtie2_species ${merge_midas_outdir}/bt2_indexes/pangenomes.species \
     ${midas_outdir} &> ${logs_dir}/xx_genes_${num_cores}_w_bowtie2.log"


echo "Testing Across-Samples CNV Module"
midas merge_genes --samples_list ${pool_fp} --midasdb_name ${midas_dbname} --midasdb_dir ${midas_db} \
     --num_cores ${num_cores} --sample_counts 2 ${merge_midas_outdir} \
     --cluster_level_in 99 --genome_depth 0.4 \
     &> ${logs_dir}/merge_genes_${num_cores}.log


echo "MIDASv3 Unit Testing SUCCESS"
