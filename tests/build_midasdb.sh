#! /usr/bin/bash
set -e

basedir=`pwd`
testdir="${basedir}/tests"
echo ${testdir}

db_dir="$testdir/$db_dir/local_testdb"

echo "Build MIDAS-DB locally $db_dir"

python3 -m iggtools annotate_genome --species all --midasdb_name testdb --midasdb_dir $db_dir  --debug --force
python3 -m iggtools build_midasdb --generate_gene_feature --genomes all --midasdb_name testdb --midasdb_dir $db_dir --debug --force


python3 -m iggtools infer_markers --genomes all --midasdb_name testdb --midasdb_dir $db_dir --debug --force
python3 -m iggtools build_midasdb --build_markerdb --midasdb_name testdb --midasdb_dir $db_dir --debug --force


python3 -m iggtools build_pangenome --species all --midasdb_name testdb --midasdb_dir $db_dir --debug --force
python3 -m iggtools build_midasdb --generate_cluster_info --species all --midasdb_name testdb --midasdb_dir $db_dir --debug --force


myArray=(100000 500000 1000000)
for chunksize in ${myArray[@]}; do
  python3 -m iggtools compute_chunks --chunk_type genes --chunk_size $chunksize --species all --midasdb_name testdb --midasdb_dir $db_dir --debug --force
  python3 -m iggtools compute_chunks --chunk_type run_snps --chunk_size $chunksize --species all --midasdb_name testdb --midasdb_dir $db_dir --debug --force
done

python3 -m iggtools compute_chunks --chunk_type merge_snps --chunk_size 500000 --species all --midasdb_name testdb --midasdb_dir $db_dir --debug --force

python3 -m iggtools build_bowtie2db --midasdb_name testdb --midasdb_dir $db_dir --species_list 117086,117088 --bt2_indexes_dir $db_dir/bt2_indexes --bt2_indexes_name repgenomes
python3 -m iggtools build_bowtie2db --midasdb_name testdb --midasdb_dir $db_dir --species_list 117086,117088 --bt2_indexes_dir $db_dir/bt2_indexes --bt2_indexes_name pangenomes
