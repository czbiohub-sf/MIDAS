#!/bin/bash                         #-- what is the language of this shell
#                                   #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o /wynton/group/pollard/czhao/midasdb_wis/pipeline_out                       #-- output directory (fill in)
#$ -e /wynton/group/pollard/czhao/midasdb_wis/pipeline_err                       #-- error directory (fill in)
#$ -V                               #-- pass all environment variables
#$ -cwd                             #-- tell the job that it should start in your working directory
##$ -r y                            #-- tell the system that if a job crashes, it should be restarted
##$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -pe smp 8
#$ -l mem_free=5G                   #-- submits on nodes with enough free memory (required)
#$ -l arch=lx-amd64                 #-- SGE resources (CPU type)
#$ -l scratch=20G                   #-- SGE resources (home and scratch disks)
#$ -l h_rt=00:30:00                #-- runtime limit (see above; this requests 24 hours)
#$ -t 1-2                        #-- Array job: submit XX jobs to the job queues at once.
#$ -tc 2                         #-- specify the maximum number of concurrent tasks


module load Sali anaconda
conda activate midas2v1.07

echo "Hello world, I’m running on node ${HOSTNAME}"
JB_LAB=batch_${SGE_TASK_ID}
TREADS=${NSLOTS:-1}

GNUTIME="/wynton/home/pollard/czhao/local/bin/time-1.9/bin/time"
#MIDASDB_DIR="/pollard/data/midas2-db/midas2db-wis-2023"
MIDASDB_DIR="/pollard/data/midas2-db/toy_wis"

CWDIR="/wynton/home/pollard/czhao/midasdb_wis"
#GB_IN="${CWDIR}/sge_jobs/step1_annotate/${JB_LAB}"
GB_IN="${CWDIR}/sge_jobs/${JB_LAB}"

#### global scratch for outputs
GB_OUT="/wynton/scratch/czhao/midasdb_wis/step1_annotate/${JB_LAB}"
rm -rf ${GB_OUT}
mkdir -p $GB_OUT

while IFS= read -r genome_id
do
  GENOME_DIR="${GB_OUT}/${genome_id}"
  TEMP_DIR=$SPECIES_DIR/temp
  if [ ! -d $GENOME_DIR/temp ]; then
    mkdir -p $GENOME_DIR/temp
  fi

  ${GNUTIME} -v midas2 annotate_genome --genomes ${genome_id} \
    --midasdb_name newdb --midasdb_dir ${MIDASDB_DIR} \
    --scratch_dir ${GENOME_DIR} --num_threads ${TREADS} --debug &> ${TEMP_DIR}/annotate_time.log

  ${GNUTIME} -v midas2 infer_markers --genomes ${genome_id} \
    --midasdb_name newdb --midasdb_dir ${MIDASDB_DIR} \
    --debug --num_threads ${TREADS} &> ${TEMP_DIR}/infer_time.log

  ${GNUTIME} -v midas2 build_midasdb --generate_gene_feature --genomes ${genome_id} \
    --midasdb_name newdb --midasdb_dir ${MIDASDB_DIR} \
    --debug --num_threads ${TREADS} &> ${TEMP_DIR}/gene_feature_time.log

done < $GB_IN



SUB_DIR=${CWDIR}
if ls -A ${SUB_DIR}/ | grep core; then
  rm -v ${SUB_DIR}/core.*
else
  echo “no core dumping detected”
fi

conda deactivate

## End-of-job summary, if running as a job
#[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"          # This is useful for debugging and usage purposes
