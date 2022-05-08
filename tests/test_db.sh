#! /usr/bin/bash
set -e

basedir=`pwd`
testdir="${basedir}/tests/midasdb_newdb"
echo ${testdir}

db_name="testdb"
db_dir="$testdir"

echo "START MIDAS 2.0 Database Build Testing"
echo "testdb is built locally at $db_dir"

#rm -rf $db_dir

echo "Annotate Genomes and Repgenome"
midas2 annotate_genome --species all --midasdb_name $db_name --midasdb_dir $db_dir --debug --force
midas2 build_midasdb --generate_gene_feature --genomes all --midasdb_name $db_name --midasdb_dir $db_dir --debug --force

echo "Infer SGC genes and Build marker DB"
midas2 infer_markers --genomes all --midasdb_name $db_name --midasdb_dir $db_dir --debug --force
midas2 build_midasdb --build_markerdb --midasdb_name $db_name --midasdb_dir $db_dir --debug --force

echo "Build Pangenomes"
midas2 build_pangenome --species all --midasdb_name $db_name --midasdb_dir $db_dir --debug --force
midas2 build_midasdb --generate_cluster_info --species all --midasdb_name $db_name --midasdb_dir $db_dir --debug --force


echo "Compute Chunks"
midas2 compute_chunks --chunk_type genes --chunk_size 50000 --species all --midasdb_name $db_name --midasdb_dir $db_dir --debug --force
midas2 compute_chunks --chunk_type run_snps --chunk_size 1000000 --species all --midasdb_name $db_name --midasdb_dir $db_dir --debug --force
midas2 compute_chunks --chunk_type merge_snps --chunk_size 500000 --species all --midasdb_name $db_name --midasdb_dir $db_dir --debug --force

echo "SUCCESS FINISH MIDAS 2.0 Database Build"
