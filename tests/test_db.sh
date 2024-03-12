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

db_name="newdb"
db_dir="$testdir/midasdb_newdb"

echo "START MIDAS2 Database Build Testing"
echo "MIDASDB $db_name is built locally at $db_dir"

rm -rf $db_dir

bash /pollard/home/czhao/projects/MIDAS2/bin/run_mge.sh /pollard/home/czhao/projects/MIDAS2/tests/midasdb_newdb 4 117086echo "Make copy of custom collection of genomes"
cp -r ${testdir}/mags $db_dir

scratch_dir="."

echo "Annotate Genomes and Repgenome"
midas2 annotate_genome --species all --midasdb_name $db_name --midasdb_dir $db_dir --debug --force --scratch_dir ${scratch_dir} -t ${num_cores}
midas2 build_midasdb --generate_gene_feature --genomes all --midasdb_name $db_name --midasdb_dir $db_dir --debug --force --scratch_dir ${scratch_dir} -t ${num_cores}

echo "Infer SGC genes and Build marker DB"
midas2 infer_markers --genomes all --midasdb_name $db_name --midasdb_dir $db_dir --debug --force -t ${num_cores}
midas2 build_midasdb --build_markerdb --midasdb_name $db_name --midasdb_dir $db_dir --debug --force --scratch_dir ${scratch_dir} -t ${num_cores}

echo "Build Pangenome"
midas2 build_pangenome --species all --midasdb_name $db_name --midasdb_dir $db_dir --debug --force --scratch_dir ${scratch_dir} -t ${num_cores} --recluster
bash $basedir/bin/pipeline.sh 117086 $db_dir/pangenomes 2 100000 $basedir/bin
bash $basedir/bin/pipeline.sh 117088 $db_dir/pangenomes 2 100000 $basedir/bin
midas2 recluster_centroids --species all --midasdb_name $db_name --midasdb_dir $db_dir --debug --force --scratch_dir ${scratch_dir} -t ${num_cores}
midas2 augment_pangenome --species all --midasdb_name $db_name --midasdb_dir $db_dir --debug --force --scratch_dir ${scratch_dir} -t ${num_cores}

echo "NEED TO INSTALL GENOMAD, EGGNOG MAPPER, RESFINDER AND MEFINDER"
grep -v genome $db_dir/genomes.tsv | cut -f1,2 | xargs -Ixx -l -P 1 bash -c 'echo bash '"$basedir/bin/run_genomad.sh $midasdb_dir $num_cores"' $0 $1'
grep -v genome $db_dir/genomes.tsv | cut -f1,2 | xargs -Ixx -l -P 1 bash -c 'echo bash '"$basedir/bin/run_mge.sh $midasdb_dir $num_cores"' $0 $1'
grep -v genome $db_dir/genomes.tsv | cut -f2 | sort | uniq | xargs -Ixx -l -P 1 bash -c 'echo bash '"$basedir/bin/run_eggnog.sh $midasdb_dir $num_cores"' $0'
cp -r $testdir/pangenomes_annotation $db_dir/pangenomes_annotation


echo "Enhance Pangenome"
midas2 annotate_pangenome --species all --midasdb_name $db_name --midasdb_dir $db_dir --debug --force --scratch_dir ${scratch_dir} -t ${num_cores}
midas2 enhance_pangenome --species all --midasdb_name $db_name --midasdb_dir $db_dir --debug --force --scratch_dir ${scratch_dir} -t ${num_cores}

echo "Compute Chunks"
midas2 compute_chunks --chunk_type run_snps --chunk_size 1000000 --species all --midasdb_name $db_name --midasdb_dir $db_dir --debug --force -t ${num_cores}
midas2 compute_chunks --chunk_type merge_snps --chunk_size 500000 --species all --midasdb_name $db_name --midasdb_dir $db_dir --debug --force -t ${num_cores}

echo "SUCCESS FINISH MIDAS2 Database Build"
