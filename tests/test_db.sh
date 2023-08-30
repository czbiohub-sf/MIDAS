#! /usr/bin/bash
set -e

basedir=`pwd`
testdir="${basedir}/tests"
echo ${testdir}

db_name="newdb"
db_dir="$testdir/midasdb_newdb"

echo "START MIDAS2 Database Build Testing"
echo "MIDASDB $db_name is built locally at $db_dir"

rm -rf $db_dir

echo "Make copy of custom collection of genomes"
cp -r ${testdir}/mags $db_dir

scratch_dir="."

echo "Annotate Genomes and Repgenome"
midas2 annotate_genome --species all --midasdb_name $db_name --midasdb_dir $db_dir --debug --force --scratch_dir ${scratch_dir}
midas2 build_midasdb --generate_gene_feature --genomes all --midasdb_name $db_name --midasdb_dir $db_dir --debug --force --scratch_dir ${scratch_dir}

echo "Infer SGC genes and Build marker DB"
midas2 infer_markers --genomes all --midasdb_name $db_name --midasdb_dir $db_dir --debug --force
midas2 build_midasdb --build_markerdb --midasdb_name $db_name --midasdb_dir $db_dir --debug --force --scratch_dir ${scratch_dir}

echo "Build Pangenomes"
midas2 build_pangenome --species all --midasdb_name $db_name --midasdb_dir $db_dir --debug --force --scratch_dir ${scratch_dir}
midas2 build_midasdb --generate_cluster_info --species all --midasdb_name $db_name --midasdb_dir $db_dir --debug --force --scratch_dir ${scratch_dir}


echo "Compute Chunks"
midas2 compute_chunks --chunk_type genes --chunk_size 50000 --species all --midasdb_name $db_name --midasdb_dir $db_dir --debug --force
midas2 compute_chunks --chunk_type run_snps --chunk_size 1000000 --species all --midasdb_name $db_name --midasdb_dir $db_dir --debug --force
midas2 compute_chunks --chunk_type merge_snps --chunk_size 500000 --species all --midasdb_name $db_name --midasdb_dir $db_dir --debug --force

echo "SUCCESS FINISH MIDAS2 Database Build"
