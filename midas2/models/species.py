#!/usr/bin/env python3
import os
import json
from math import floor
from collections import defaultdict
from operator import itemgetter

from midas2.common.utils import InputStream, OutputStream, command, select_from_tsv
from midas2.common.utilities import scan_fasta, scan_cluster_info
from midas2.params.schemas import fetch_cluster_xx_info_schema


class Species:
    """ Base class for species """
    def __init__(self, species_id):
        self.id = species_id

        # SNPs chunk
        self.contigs_fp = None
        self.contigs = defaultdict(dict)
        self.chunks_of_sites_fp = None
        self.num_of_snps_chunks = None
        self.max_contig_length = None

        # Genes
        self.clusters_info_fp = {} # Initialize an empty dictionary for cluster_xx_info
        self.pangenome_size = {}
        self.clusters_info = {} # Initialize an empty dictionary for cluster_xx_info
        self.list_of_markers = {}
        self.clusters_map = {}

        # MERGE Flow: select species
        self.list_of_samples = [] # relevant samples for given species
        self.samples_count = 0
        self.list_of_samples_depth = [] # mean genome coverage

        # Merge SNPs
        self.gene_feature_fp = None
        self.gene_seq_fp = None


    def set_clusters_info_fp(self, midas_db, xx):
        # Dynamically create or update a cluster attribute
        # The attribute name is constructed dynamically based on xx
        species_id = self.id
        self.clusters_info_fp[xx] = midas_db.get_target_layout("pangenome_cluster_xx", False, species_id, "", xx)


    def get_clusters_info_fp(self, xx):
        # Dynamically retrieve the cluster attribute
        return self.clusters_info_fp.get(xx, None)


    def get_clusters_info(self, xx):
        # Dynamically retrieve the cluster attribute
        if xx == "99":
            cxx_info = scan_cluster_info(self.get_clusters_info_fp(xx), xx)
        else:
            cxx_info = scan_cluster_info(self.get_clusters_info_fp(xx), xx, fetch_cluster_xx_info_schema(xx))

        list_of_markers = list(set(r[f'centroid_{xx}_marker_id'] for r in cxx_info.values() if r[f'centroid_{xx}_marker_id']))

        self.clusters_info[xx] = cxx_info
        self.pangenome_size[xx] = len(cxx_info)
        self.list_of_markers[xx] = list_of_markers


    def get_cluster_map(self, xx_in, xx_out):
        c99_info = self.clusters_info['99']
        cid_map = {}
        for row in c99_info.values():
            cxx_in = f"centroid_{xx_in}"
            cxx_out = f"centroid_{xx_out}"
            if cxx_in not in cid_map:
                cid_map[row[cxx_in]] = row[cxx_out]
        self.clusters_map[xx_in] = cid_map


    def compute_snps_chunks(self, midas_db, chunk_size, workflow):
        """ The structure of the chunks depends on the repgenome and chunk size """

        species_id = self.id
        genome_id = midas_db.uhgg.fetch_repgenome_id(species_id)

        if workflow == "run":
            chunk_filename = "chunks_sites_run"
        if workflow == "merge":
            chunk_filename = "chunks_sites_merge"

        local_file = midas_db.get_target_layout(chunk_filename, False, species_id, genome_id, chunk_size)
        contigs_fp = midas_db.get_target_layout("representative_genome", False, species_id, genome_id)

        if not os.path.exists(local_file):
            if workflow == "run":
                chunks_of_sites = design_run_snps_chunks(species_id, contigs_fp, chunk_size)
            if workflow == "merge":
                chunks_of_sites = design_merge_snps_chunks(species_id, contigs_fp, chunk_size)
            write_chunks_cache(chunks_of_sites, local_file)

        chunks_of_sites = load_chunks_cache(local_file)
        _, _, number_of_chunks, max_contig_length = chunks_of_sites[-1]

        self.chunks_of_sites_fp = local_file
        self.num_of_snps_chunks = number_of_chunks
        self.max_contig_length = max_contig_length

        if workflow == "run":
            self.contigs_fp = contigs_fp

        if workflow == "merge":
            chunks_contigs = {}
            for chunk_list in chunks_of_sites.values():
                _, chunk_id, contig_id, list_of_contigs = chunk_list[0][:4]
                loc_fp = midas_db.get_target_layout("chunks_contig_lists", False, species_id, chunk_id, chunk_size)

                if contig_id == -1:
                    chunks_contigs[chunk_id] = loc_fp
                    if not os.path.exists(loc_fp):
                        command(f"mkdir -p {os.path.dirname(loc_fp)}")
                        with OutputStream(loc_fp) as stream:
                            stream.write("\n".join(list_of_contigs) + "\n")

            self.chunks_contigs = chunks_contigs
            self.gene_feature_fp = midas_db.fetch_file("annotation_genes", species_id)
            self.gene_seq_fp = midas_db.fetch_file("annotation_ffn", species_id)

        return chunks_of_sites


    def get_repgenome(self, midas_db):
        species_id = self.id
        genome_id = midas_db.get_repgenome_id(species_id)
        self.contigs = scan_fasta(midas_db.get_target_layout("representative_genome", False, species_id, genome_id))


    def fetch_contigs_ids(self):
        list_of_contig_ids = []
        with InputStream(self.contigs_fp, "grep \'>\'") as stream:
            for line in stream:
                list_of_contig_ids.append(line.strip("\n")[1:])
        stream.ignore_errors()
        return list_of_contig_ids


    def fetch_samples_names(self):
        return [sample.sample_name for sample in self.list_of_samples]


    def fetch_samples_depth(self):
        self.list_of_samples_depth = [sample.profile[self.id]["mean_depth"] for sample in self.list_of_samples]


def parse_species(args):
    species_list = []
    if args.species_list:
        if os.path.exists(args.species_list):
            with InputStream(args.species_list) as stream:
                for line in stream:
                    species_list.append(line.strip())
        else:
            species_list = args.species_list.split(",")
    return species_list


def filter_species(profile_fp, select_by, select_threshold, species_list=None):
    select_by = select_by.split(",")
    select_threshold = select_threshold.split(",")
    nargs = len(select_by)
    assert len(select_by) == len(select_threshold)

    dict_of_species = {}
    if "median_marker_coverage" in select_by:
        column_names = list(set(["species_id", "median_marker_coverage"] + select_by))
    else:
        column_names = ["species_id"] + select_by

    with InputStream(profile_fp) as stream:
        for record in select_from_tsv(stream, selected_columns=column_names, result_structure=dict):
            if species_list and record["species_id"] not in species_list:
                continue
            if sum([1 if float(record[select_by[i]]) > float(select_threshold[i]) else 0 for i in range(nargs)]) == nargs:
                dict_of_species[record["species_id"]] = float(record[column_names[1]])

    # Sort species in descending order of first column name (median_marker_coverage)
    species_ids = [k for k, v in sorted(dict_of_species.items(), key=itemgetter(1), reverse=True)]
    return species_ids


def load_chunks_cache(chunk_cache_fp):
    assert os.path.exists(chunk_cache_fp), f"{chunk_cache_fp} doesn't exit"
    with InputStream(chunk_cache_fp) as stream:
        chunks_dict = json.load(stream)
        chunks_dict = {int(k):v for k, v in chunks_dict.items()} # conver back to int key
    return chunks_dict


def write_chunks_cache(chunks_of_sites, local_file):
    command(f"mkdir -p {os.path.dirname(local_file)}")
    with OutputStream(local_file) as stream:
        json.dump(chunks_of_sites, stream)


def partition_contigs_into_chunks(unassigned_contigs, chunk_size, chunk_id):
    """ Partition short, unassigned contigs into subsets/chunks.
        Similar to the problem of partition to K equal sum subsets """

    # Sort contigs by descending order of contig length
    sorted_contigs = {cid:cc["contig_length"] for cid, cc in sorted(unassigned_contigs.items(), key=lambda x: x[1]["contig_length"], reverse=True)}
    list_of_contigs_id = list(sorted_contigs.keys())
    list_of_contigs_length = list(sorted_contigs.values())

    # Initialize two pointers
    istart = 0
    jstart = len(list_of_contigs_length)-1
    prev_jstart = jstart + 1
    subset_of_contigs = defaultdict(dict)
    list_of_partitoned_contigs = [] # for valication purpose

    while istart <= jstart:
        # Add one long contig to current chunk
        curr_chunk_length = list_of_contigs_length[istart]
        # Iteratively add smaller contigs until exceed chunk size
        while curr_chunk_length + sum(list_of_contigs_length[jstart:prev_jstart]) <= chunk_size and istart < jstart:
            jstart = jstart - 1

        # Collect all the added shorted contigs
        added_clens = list_of_contigs_length[jstart+1:prev_jstart]
        curr_clens = [curr_chunk_length] + added_clens
        curr_chunk_length += sum(added_clens)

        # Record the list of contig_ids assigned to current chunk_id
        curr_cids = [list_of_contigs_id[istart]] + list_of_contigs_id[jstart+1:prev_jstart]
        subset_of_contigs[chunk_id] = {
            "chunk_id": chunk_id,
            "contigs_id": curr_cids,
            "chunk_length": curr_chunk_length,
            "list_of_contigs_length": curr_clens
        }

        # Update the pointer and chunk_id
        list_of_partitoned_contigs = list_of_partitoned_contigs + curr_clens
        istart = istart + 1
        prev_jstart = jstart + 1
        chunk_id += 1

    assert len(list_of_partitoned_contigs) == len(list_of_contigs_length)
    assert set(list_of_partitoned_contigs) == set(list_of_contigs_length)

    return (subset_of_contigs, chunk_id)


def design_run_snps_chunks(species_id, contigs_file, chunk_size):
    """ Given the Genome and chunk_size, the structure of the chunks are the same.
        Each chunk is indexed by (species_id, chunk_id) """

    # The structure of the chunks depends on the representative genome sequences
    contigs = scan_fasta(contigs_file)

    # Start with full chunks
    chunk_id = 0
    chunks_of_sites = defaultdict(list) # list of tuples
    unassigned_contigs = defaultdict(dict)
    max_contig_length = 0

    for contig in contigs.values():
        contig_length = contig["length"]
        if contig_length > max_contig_length:
            max_contig_length = contig_length

        contig_id = contig["id"]
        # left closed, right open
        if contig_length < chunk_size:
            unassigned_contigs[contig_id] = {"contig_id": contig_id,
                                             "contig_start": 0,
                                             "contig_end": contig_length,
                                             "contig_length": contig_length,
                                             "compute_reads": True}
        else:
            number_of_full_chunks = floor(contig_length/chunk_size)
            for ni, ci in enumerate(range(0, contig_length, chunk_size)):
                if ni == number_of_full_chunks: # last chunk
                    unassigned_contigs[contig_id] = {"contig_id": contig_id,
                                                     "contig_start": ci,
                                                     "contig_end": contig_length,
                                                     "contig_length": contig_length - ci,
                                                     "compute_reads": False}
                else:
                    count_flag = ni == 0 # first chunk
                    chunks_of_sites[chunk_id] = [(species_id, chunk_id, contig_id, ci, ci+chunk_size, count_flag, 0)]
                    chunk_id += 1

    if unassigned_contigs:
        # Partition unassigned short contigs into subsets
        subset_of_contigs, chunk_id = partition_contigs_into_chunks(unassigned_contigs, chunk_size, chunk_id)

        # Add the partitioned subsets to chunks
        for chunk_dict in subset_of_contigs.values():
            _chunk_id = chunk_dict["chunk_id"]
            list_of_contigs = chunk_dict["contigs_id"]
            for _cidx, _cid in enumerate(list_of_contigs):
                cstart = unassigned_contigs[_cid]["contig_start"]
                cend = unassigned_contigs[_cid]["contig_end"]
                cflag = unassigned_contigs[_cid]["compute_reads"]
                chunks_of_sites[_chunk_id].append((species_id, _chunk_id, _cid, cstart, cend, cflag, _cidx))
        assert chunk_id == _chunk_id+1

    # Finally the merge jobs
    number_of_chunks = chunk_id
    chunks_of_sites[-1] = (species_id, -1, number_of_chunks, max_contig_length)

    return chunks_of_sites


def design_merge_snps_chunks(species_id, contigs_file, chunk_size):

    contigs = scan_fasta(contigs_file)

    # Start with full chunks
    chunk_id = 0
    chunks_of_sites = defaultdict(list)
    unassigned_contigs = defaultdict(dict)
    max_contig_length = 0

    for contig in contigs.values():
        contig_length = contig["length"]
        if contig_length > max_contig_length:
            max_contig_length = contig_length

        contig_id = contig["id"]
        # left closed, right open
        if contig_length < chunk_size:
            unassigned_contigs[contig_id] = {"contig_id": contig_id,
                                             "contig_start": 0,
                                             "contig_end": contig_length,
                                             "contig_length": contig_length}
        else:
            number_of_full_chunks = floor(contig_length/chunk_size)
            for ni, ci in enumerate(range(0, contig_length, chunk_size)):
                if ni == number_of_full_chunks - 1: # last full chunk: carry over
                    chunks_of_sites[chunk_id] = [(species_id, chunk_id, contig_id, ci, contig_length)]
                    chunk_id += 1
                    break
                else:
                    chunks_of_sites[chunk_id] = [(species_id, chunk_id, contig_id, ci, ci+chunk_size)]
                    chunk_id += 1

    if unassigned_contigs:
        # Partition unassigned short contigs into subsets
        dict_of_chunks, chunk_id = partition_contigs_into_chunks(unassigned_contigs, chunk_size, chunk_id)
        # Add the partitioned subsets to chunks
        for chunk_dict in dict_of_chunks.values():
            _chunk_id = chunk_dict["chunk_id"]
            chunks_of_sites[_chunk_id].append((species_id, _chunk_id, -1, chunk_dict["contigs_id"]))
        assert chunk_id == _chunk_id+1

    number_of_chunks = chunk_id
    chunks_of_sites[-1] = (species_id, -1, number_of_chunks, max_contig_length)

    return chunks_of_sites
