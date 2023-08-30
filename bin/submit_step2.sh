#!/bin/bash                         #-- what is the language of this shell
#                                   #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o /wynton/scratch/czhao/midasdb_wis/step2_out                       #-- output directory (fill in)
#$ -e /wynton/scratch/czhao/midasdb_wis/step2_err                       #-- error directory (fill in)
#$ -V                               #-- pass all environment variables
#$ -cwd                             #-- tell the job that it should start in your working directory
##$ -r y                            #-- tell the system that if a job crashes, it should be restarted
##$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -pe smp 8
#$ -l mem_free=8G                   #-- submits on nodes with enough free memory (required)
#$ -l arch=lx-amd64                 #-- SGE resources (CPU type)
#$ -l scratch=40G                   #-- SGE resources (home and scratch disks)
#$ -l h_rt=24:00:00                #-- runtime limit (see above; this requests 24 hours)
#$ -t 1-2                        #-- Array job: submit XX jobs to the job queues at once.
#$ -tc 2                         #-- specify the maximum number of concurrent tasks


module load Sali anaconda
conda activate midas2v1.07

echo "Hello world, I’m running on node ${HOSTNAME}"
JB_LAB=species_${SGE_TASK_ID}
TREADS=${NSLOTS:-1}

GNUTIME="/wynton/home/pollard/czhao/local/bin/time-1.9/bin/time"
MIDASDB_DIR="/wynton/scratch/czhao/midas2db-wis-2023"

CWDIR="/wynton/home/pollard/czhao/midasdb_wis"
SCRIPTDIR="${CWDIR}/MIDAS2/bin"
GB_IN="${CWDIR}/sge_jobs/step2_pangenome/small/${JB_LAB}"

# local scratch disk as output directory
OUTDIR="/scratch/czhao/step2_pangenome/${JB_LAB}"
if [ ! -d $OUTDIR ]; then
  mkdir -p $OUTDIR
fi

cd "$OUTDIR" # Use a temporary working directory

LOGDIR="/wynton/scratch/czhao/midasdb_wis/step2_log/${JB_LAB}"
mkdir -p ${LOGDIR}

SCRATCH_DIR="$OUTDIR/scratch"
if [ ! -d $SCRATCH_DIR ]; then
  mkdir -p $SCRATCH_DIR
fi

while IFS= read -r species_id
do
  ${GNUTIME} -v python -m midas2 build_pangenome \
    --midasdb_name newdb --midasdb_dir ${MIDASDB_DIR} \
    --species ${species_id} --num_threads ${TREADS} \
    --scratch_dir ${SCRATCH_DIR} \
    --debug --recluster &> ${LOGDIR}/${species_id}/pangenome_time.log

  ${GNUTIME} -v bash ${SCRIPTDIR}/pipeline.sh ${species_id} \
      ${SPECIES_DIR} ${TREADS} 64000 \
    &> ${LOGDIR}/${species_id}/pipeline_time.log

  ${GNUTIME} -v python -m midas2 recluster_centroids \
    --midasdb_name newdb --midasdb_dir ${MIDASDB_DIR} \
    --species ${species_id} --num_threads ${TREADS} \
    --scratch_dir ${SCRATCH_DIR} \
    --debug &> ${LOGDIR}/${species_id}/recluster_time.log

  ${GNUTIME} -v python -m midas2 build_midasdb --generate_cluster_info \
    --midasdb_name newdb --midasdb_dir ${MIDASDB_DIR} \
    --species ${species_id} --num_threads ${TREADS} \
    --scratch_dir ${SCRATCH_DIR} \
    --debug &> ${LOGDIR}/${species_id}/cluster_info_time.log

  cp -r ${SCRATCH_DIR}/buildpan/${species_id}/temp/cdhit ${MIDASDB_DIR}/pangenomes/${species_id}/temp
  cp ${LOGDIR}/${species_id}/*.log ${MIDASDB_DIR}/pangenomes/${species_id}/temp

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
