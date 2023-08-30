# Chunyu Zhao 2023-06-14
# Dependency:
# seqkit
# awk

#! /usr/bin/bash

if [ $# -ne 4 ]; then
    echo "Usage: $0 SPECIES SPECIEDIR THREADS MEM"
    exit 1
fi

#######
script_dir="/wynton/home/pollard/czhao/midasdb_wis/MIDAS2/bin" #<----

species_id="$1"
species_dir="$2"
total_threads="$3"
total_mem="$4" #500000

vsearch_threads=8
vsearch_jobs=$((total_threads / vsearch_threads))
thread2=$((total_threads / 2))

####### INPUTS: global scratch directory
vsearch_dir="${species_dir}/temp/vsearch"
genes_ffn="${vsearch_dir}/genes.ffn"
genes_len="${vsearch_dir}/genes.len"
genes_info="${vsearch_dir}/gene_info.txt"
centroids_clean_ffn="${vsearch_dir}/centroids.99.clean.ffn"
centroids_ambig_ffn="${vsearch_dir}/centroids.99.ambiguous.ffn"


####### OUTPUTS
out_dir="${species_dir}/temp/cdhit"
mkdir -p ${out_dir}

members_dir="${out_dir}/step1_members"
vsearch_dir="$out_dir/step2_vsearch"
cdhit_dir="$out_dir/step3_cdhit"
info_dir="$out_dir/step4_info"
mkdir -p $members_dir
mkdir -p $vsearch_dir
mkdir -p ${cdhit_dir}
mkdir -p $info_dir


####### Replace centroid_99 with ambiguous with new centroids
gene_info_vsearch="${info_dir}/gene_info_vsearch.tsv"
vsearch_centroids_ffn="${info_dir}/vsearch_centroids.ffn"
if [[ -s ${centroids_ambig_ffn} ]]; then
  if [[ ! -e "${gene_info_vsearch}" ]]; then
    # Get the list fo centroids99 with ambiguous bases:
    ambiguous_list="${members_dir}/list_of_ambiguous_centroids"
    grep ">" ${centroids_ambig_ffn} | awk 'sub(/>/, "")' > $ambiguous_list

    # For each ambiguous centroids, get all the member sequences.
    # Two intermediate files: members.id and members.ffn.
    cat ${ambiguous_list} | \
      xargs -Ixx -P ${thread2} bash -c "bash ${script_dir}/get_members.sh xx $genes_ffn $genes_info $members_dir/xx.mems.ffn"

    # Get list of centroids with more than one members
    multimember_centroids="${members_dir}/centroids_with_multimembers"
    find ${members_dir} -name '*mems.ffn' -size +0c | awk 'sub(/\.mems\.ffn/, "")' | awk -F/ '{print $NF}' > ${multimember_centroids}

    # Select new centroids for the previously ambiguous centroid99 clusters
    cat ${multimember_centroids} | \
      xargs -Ixx -P ${vsearch_jobs} bash -c "vsearch --cluster_fast $members_dir/xx.mems.ffn --threads ${vsearch_threads} --quiet --id 0.99 --centroids $vsearch_dir/xx.centroids -uc $vsearch_dir/xx.clusters"

    # gather gene info for newly clustered centroids_99 via research
    cat ${multimember_centroids} | \
      xargs -Ixx -P ${thread2} bash -c "bash ${script_dir}/gather_geneinfo.sh $vsearch_dir/xx.clusters $vsearch_dir/xx.geneinfo"

    # collect related gene info changes
    gene_info_add="$info_dir/vsearch_gene_info_add.tsv"
    cat ${vsearch_dir}/*.geneinfo | tr ' ' '\t' > $gene_info_add
    cut -f2 $gene_info_add  | sort | uniq > $vsearch_dir/list_of_vsearch_centroids
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$0];next} !($2 in a)' ${ambiguous_list} ${genes_info} > ${vsearch_dir}/vsearch_gene_info_keep.tsv

    cat ${vsearch_dir}/vsearch_gene_info_keep.tsv $gene_info_add > ${gene_info_vsearch}
  fi

  if [[ ! -e "${vsearch_centroids_ffn}" ]]; then
    # Add new centroids_99 to entroids.99.clean.ffn
    seqkit grep -w 0 -f $vsearch_dir/list_of_vsearch_centroids ${genes_ffn} > "${vsearch_dir}/centroids_add.ffn"
    cat ${centroids_clean_ffn} "${vsearch_dir}/centroids_add.ffn" > ${vsearch_centroids_ffn}
  fi
else
  echo "No ambiguous centroids_99"
  cp ${genes_info} ${gene_info_vsearch}
  cp ${centroids_clean_ffn} ${vsearch_centroids_ffn}
fi


cdhit_centroids_ffn="${info_dir}/cdhit_centroids.ffn"
cdhit_centroids_tsv="${info_dir}/cdhit_centroids.tsv"
if [[ ! -e "${cdhit_centroids_tsv}" ]]; then
  ##### Input: vsearch_centroids_ffn
  pushd ${cdhit_dir}
  cp ${vsearch_centroids_ffn} vsearch_centroids.ffn

  # rename the input fasta with a short and save the mapping file
  awk '/^>/{$0=">g"++i; }1' vsearch_centroids.ffn > vsearch_centroids_renamed.ffn
  paste -d '\t' <(grep "^>" vsearch_centroids.ffn  | cut -c 2-) <(grep "^>" vsearch_centroids_renamed.ffn | cut -c 2-) > mapping.txt

  # detect partial gene duplication + reverse complement
  cd-hit-est -i vsearch_centroids_renamed.ffn -c 1 -T ${total_threads} -aS 0.9 -G 0 -g 1 -AS 180 -M ${total_mem} -o cdhit_centroids.ffn

  # make sure all sequences are in upper letter
  awk '{ if ($0 ~ /^>/) {print $0} else {print toupper($0)}}' cdhit_centroids.ffn > cdhit_centroids_upper.ffn
  # get vsearch centroids_99 to cdhit centorids_99 mapping
  bash ${script_dir}/parse_cdhit_cluster.sh cdhit_centroids.ffn.clstr cdhit_centroids.ffn.clstr.tsv
  # 2023-07-24 cd-hit clstr run of out space error: solution: rename the original input sequences

  # get cdhit generated centroids
  awk -v OFS='\t' '$2 == 0 {print $1, $4, $4}' cdhit_centroids.ffn.clstr.tsv > new_centroids.tsv
  # get cdhit generated members
  awk -v OFS='\t' '$2 != 0 {print $1, $4}' cdhit_centroids.ffn.clstr.tsv > new_members.tsv
  # add corresponding centroids to members
  join -t $'\t' <(sort -k1,1 new_members.tsv) <(sort -k1,1 new_centroids.tsv | cut -f1-2) > new_members_w_centroids.tsv

  cat new_centroids.tsv new_members_w_centroids.tsv | sort -k1,2n > cdhit_centroids.tsv

  # 2023-06-30 CD-HIT Bug: some centroids from the FFN was mis-labled as members in the CLSTR file, while
  # shorter member was labled as controids in the CLSTR file.
  # Check if inconsistency exist.
  grep ">" cdhit_centroids_upper.ffn | cut -d'>' -f2 > list_of_cdhit_centroids_by_ffn
  cut -f2 new_centroids.tsv | sort | uniq > list_of_cdhit_centroids_by_tsv
  # If detect inconsistent centroids assignments:
  mislabled_centroids=$(grep -Fvwf list_of_cdhit_centroids_by_tsv list_of_cdhit_centroids_by_ffn)
  if [ -n "$mislabled_centroids" ]; then
    # loop over the problematic centroids
    while read -r xx; do
      echo "Mis-labled: $xx"
      bash ${script_dir}/replace_centroids.sh $xx ${cdhit_dir}/cdhit_centroids.tsv
    done <<< "$mislabled_centroids"
  fi

  ####### OUTPUT: ${cdhit_dir}/cdhit_centroids.tsv AND ${cdhit_dir}/cdhit_centroids_upper.ffn
  awk 'BEGIN{while(getline <"mapping.txt") dict[">"$2] = ">"$1} {if ($0 in dict) print dict[$0]; else print $0}' cdhit_centroids_upper.ffn > ${cdhit_centroids_ffn}
  awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$2]=$1; next} $2 in a {$2=a[$2]} 1' mapping.txt cdhit_centroids.tsv > cdhit_centroids_c2.tsv
  awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$2]=$1; next} $3 in a {$3=a[$3]} 1' mapping.txt cdhit_centroids_c2.tsv > ${cdhit_centroids_tsv}

  popd
fi


gene_info_cdhit="${info_dir}/gene_info_cdhit.tsv"
if [[ ! -e ${gene_info_cdhit} ]]; then
  # target: geneid - centroids99.cdhit
  join -t $'\t' <(awk 'BEGIN{OFS=FS="\t"}{print $2, $1}' ${gene_info_vsearch} | sort -k1,1) <(awk 'BEGIN{OFS=FS="\t"}{print $2, $3}' ${cdhit_centroids_tsv} | sort -k1,1 ) > ${info_dir}/cdhit_gene_centroid_mapping.tsv
  awk -v OFS='\t' '{print $2, $3}' ${info_dir}/cdhit_gene_centroid_mapping.tsv > ${gene_info_cdhit} #<-- this is no headers here
  cut -f1 $gene_info_cdhit > ${info_dir}/list_of_cdhit_genes
fi


is_success="${out_dir}/PIPELINE_SUCCESS"
if [[ ! -e ${is_success} ]]; then
  awk '{ if ($0 ~ /^>/) {print $0} else {print toupper($0)}}' ${cdhit_centroids_ffn} > ${out_dir}/centroids.99.ffn
  seqkit grep -w 0 -f ${info_dir}/list_of_cdhit_genes ${genes_ffn} > ${out_dir}/genes.ffn
  cp ${gene_info_cdhit} "${out_dir}/gene_info.txt"
  grep -Fwf <(cut -f1 ${gene_info_cdhit}) ${genes_len} > "${out_dir}/genes.len"

  if [ -s ${out_dir}/centroids.99.ffn ]; then
    touch $is_success
  fi
fi
