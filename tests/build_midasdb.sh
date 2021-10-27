#! /usr/bin/bash
set -e
set -x

echo "Build MIDAS-DB locally"

python -m iggtools annotate_genome --species all --midasdb_name testdb --midasdb_dir local_testdb  --debug
python -m iggtools build_midasdb --generate_gene_feature --genomes all --midasdb_name testdb --midasdb_dir local_testdb --debug


python -m iggtools infer_markers --genomes all --midasdb_name testdb --midasdb_dir local_testdb --debug
python -m iggtools build_midasdb --build_markerdb --midasdb_name testdb --midasdb_dir local_testdb --debug


python -m iggtools build_pangenome --species all --midasdb_name testdb --midasdb_dir local_testdb --debug
python -m iggtools build_midasdb --generate_cluster_info --species all --midasdb_name testdb --midasdb_dir local_testdb --debug


myArray=(50000 100000 500000 1000000 2000000)
for chunksize in ${myArray[@]}; do
  python -m iggtools compute_chunks --chunk_type genes --chunk_size $chunksize --species all --midasdb_name testdb --midasdb_dir local_testdb --debug
  python -m iggtools compute_chunks --chunk_type run_snps --chunk_size $chunksize --species all --midasdb_name testdb --midasdb_dir local_testdb --debug
done

python -m iggtools compute_chunks --chunk_type merge_snps --chunk_size 500000 --species all --midasdb_name testdb --midasdb_dir local_testdb --debug


python -m iggtools build_bowtie2db --midasdb_name testdb --midasdb_dir local_testdb --species_list 117086,117088 --bt2_indexes_dir local_testdb/bt2_indexes --bt2_indexes_name repgenomes
python -m iggtools build_bowtie2db --midasdb_name testdb --midasdb_dir local_testdb --species_list 117086,117088 --bt2_indexes_dir local_testdb/bt2_indexes --bt2_indexes_name pangenomes
