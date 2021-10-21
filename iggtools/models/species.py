#!/usr/bin/env python3
import os
import json
from math import floor
from collections import defaultdict

from iggtools.common.utils import InputStream, OutputStream, retry, command, select_from_tsv, cat_files, find_files, tsprint
from iggtools.common.utilities import scan_fasta, scan_cluster_info

from iggtools.models.midasdb import MIDAS_DB


@retry
def find_files_with_retry(f):
    return find_files(f)


class Species:
    """ Base class for species """
    def __init__(self, species_id):
        self.id = species_id

        # SNPs chunk
        self.contigs = defaultdict(dict)
        self.chunks_of_sites_fp = None
        self.num_of_snps_chunks = None
        self.max_contig_length = None

        # Genes chunk
        self.cluster_info_fp = None
        self.chunks_of_centroids_fp = None
        self.num_of_genes_chunks = None
        self.num_of_centroids = None

        # MERGE Flow: select species
        self.list_of_samples = [] # samples associate with current species
        self.samples_count = None
        self.list_of_samples_depth = [] # mean_coverage

        # Merge Genes
        self.gene_feature_file = None
        self.gene_seq_file = None


    def compute_gene_chunks(self, midas_db, chunk_size):
        """ Each chunk is indexed by <species_id, chunk_id> """

        species_id = self.id

        s3_file = midas_db.get_target_layout("chunks_centroids", True, species_id, "", chunk_size)
        local_file = midas_db.get_target_layout("chunks_centroids", False, species_id, "", chunk_size)

        if os.path.exists(local_file):
            pass
        elif find_files_with_retry(s3_file):
            midas_db.fetch_file("chunks_centroids", species_id, "", chunk_size)
        else:
            # Compute gene chunks and write to local file
            cluster_info_fp = midas_db.fetch_file("pangenome_cluster_info", species_id)
            chunks_of_centroids = design_genes_chunks(species_id, cluster_info_fp, chunk_size)
            command(f"mkdir -p {os.path.dirname(local_file)}")
            with OutputStream(local_file) as stream:
                json.dump(chunks_of_centroids, stream)

        self.chunks_of_centroids_fp = local_file
        self.get_cluster_info_fp(midas_db)

        # We load in the chunk cache only to get the num_of_gene_chunks ....
        chunks_of_centroids = load_chunks_cache(local_file)
        _, _, num_of_genes_chunks, num_of_centroids = chunks_of_centroids[-1]
        self.num_of_genes_chunks = num_of_genes_chunks
        self.num_of_centroids = num_of_centroids

        return True


    def compute_snps_chunks_run(self, midas_db, chunk_size):

        species_id = self.id
        genome_id = midas_db.uhgg.representatives[species_id]

        s3_file = midas_db.get_target_layout("chunks_sites_run", True, species_id, genome_id, chunk_size)
        local_file = midas_db.get_target_layout("chunks_sites_run", False, species_id, genome_id, chunk_size)

        if os.path.exists(local_file):
            pass
        elif find_files_with_retry(s3_file):
            midas_db.fetch_file("chunks_sites_run", species_id, "", chunk_size)
        else:
            # Compute gene chunks and write to local file
            contigs_fp = midas_db.fetch_file("representative_genome", species_id, genome_id)
            chunks_of_sites = design_run_snps_chunks(species_id, contigs_fp, chunk_size)
            command(f"mkdir -p {os.path.dirname(local_file)}")
            with OutputStream(local_file) as stream:
                json.dump(chunks_of_sites, stream)

        chunks_of_sites = load_chunks_cache(local_file)
        _, _, number_of_chunks, max_contig_length = chunks_of_sites[-1]

        self.chunks_of_sites_fp = local_file
        self.num_of_snps_chunks = number_of_chunks
        self.max_contig_length = max_contig_length

        return True


    def compute_snps_chunks_merge(self, midas_db, chunk_size):

        species_id = self.id
        genome_id = midas_db.uhgg.representatives[species_id]

        s3_file = midas_db.get_target_layout("chunks_sites_merge", True, species_id, genome_id, chunk_size)
        local_file = midas_db.get_target_layout("chunks_sites_merge", False, species_id, genome_id, chunk_size)

        if os.path.exists(local_file):
            pass
        elif find_files_with_retry(s3_file):
            midas_db.fetch_file("chunks_sites_merge", species_id, "", chunk_size)
        else:
            # Compute gene chunks and write to local file
            contigs_fp = midas_db.fetch_file("representative_genome", species_id, genome_id)
            chunks_of_sites = design_merge_snps_chunks(species_id, contigs_fp, chunk_size)
            command(f"mkdir -p {os.path.dirname(local_file)}")
            with OutputStream(local_file) as stream:
                json.dump(chunks_of_sites, stream)

        chunks_of_sites = load_chunks_cache(local_file)
        _, _, number_of_chunks, max_contig_length = chunks_of_sites[-1]

        self.chunks_of_sites_fp = local_file
        self.num_of_snps_chunks = number_of_chunks
        self.max_contig_length = max_contig_length

        return True


    def get_repgenome(self, midas_db):
        species_id = self.id
        genome_id = midas_db.get_repgenome_id(species_id)
        self.contigs = scan_fasta(midas_db.fetch_file("representative_genome", species_id, genome_id))


    def get_cluster_info_fp(self, midas_db):
        species_id = self.id
        self.cluster_info_fp = midas_db.get_target_layout("pangenome_cluster_info", False, species_id)


    def fetch_samples_names(self):
        return [sample.sample_name for sample in self.list_of_samples]


    def fetch_samples_depth(self):
        self.list_of_samples_depth = [sample.profile[self.id]["mean_coverage"] for sample in self.list_of_samples]


def load_chunks_cache(chunk_cache_fp):
    assert os.path.exists(chunk_cache_fp), f"{chunk_cache_fp} doesn't exit"
    with InputStream(chunk_cache_fp) as stream:
        chunks_dict = json.load(stream)
        chunks_dict = {int(k):v for k, v in chunks_dict.items()} # conver back to int key
    return chunks_dict


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


def filter_species(profile_fp, select_by, select_threhold, species_list=None):
    species_ids = list()
    with InputStream(profile_fp) as stream:
        for record in select_from_tsv(stream, selected_columns=["species_id", select_by], result_structure=dict):
            if species_list and record["species_id"] not in species_list:
                continue
            if float(record[select_by]) >= select_threhold: #<--
                species_ids.append(record["species_id"])
    return species_ids


def sort_list_of_species(list_of_species, rev=True):
    """ Sort list_of_species by samples_count in descending order """
    species_sorted = sorted(((sp, len(sp.list_of_samples)) for sp in list_of_species), key=lambda x: x[1], reverse=rev)
    return [sp[0] for sp in species_sorted]


def collect_chunk_pileup(sample, number_of_contigs, species_id, chunk_id, filename):
    ## Clean up sliced chunk temporary files
    headerless_sliced_path = sample.get_target_layout(filename, species_id, chunk_id)

    if number_of_contigs > 1:
        list_of_sliced_files = [sample.get_target_layout(f"{filename}_perc", species_id, chunk_id, cidx) for cidx in range(0, number_of_contigs)]
        cat_files(list_of_sliced_files, headerless_sliced_path, 20)
        for s_file in list_of_sliced_files:
            command(f"rm -rf {s_file}", quiet=True)
    else:
        sliced_file = sample.get_target_layout(f"{filename}_perc", species_id, chunk_id, 0)
        command(f"mv {sliced_file} {headerless_sliced_path}", quiet=True)


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
        subset_of_contigs[chunk_id] = {"chunk_id": chunk_id, "contigs_id": curr_cids, "chunk_length": curr_chunk_length, "list_of_contigs_length": curr_clens}
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


def design_genes_chunks(species_id, cluster_info_path, chunk_size):
    """ Each chunk is indexed by (species_id, chunk_id) """

    cluster_info = scan_cluster_info(cluster_info_path)
    num_of_centroids = len(cluster_info)

    chunk_id = 0
    genes_counter = 0
    curr_chunk_of_genes = dict()
    chunks_of_centroids = defaultdict(dict)

    for row in cluster_info.values():
        if not chunk_id*chunk_size <= genes_counter < (chunk_id+1)*chunk_size:
            chunks_of_centroids[chunk_id] = curr_chunk_of_genes
            chunk_id += 1
            curr_chunk_of_genes = defaultdict()
        curr_chunk_of_genes[row["centroid_99"]] = row["centroid_99_length"]
        genes_counter += 1

    chunks_of_centroids[chunk_id] = curr_chunk_of_genes # last chunk
    chunk_id += 1

    number_of_chunks = chunk_id
    chunks_of_centroids[-1] = (species_id, -1, number_of_chunks, num_of_centroids)

    return chunks_of_centroids


def test_chunks():
    species_id = "117086"
    chunk_size = 50000
    chunk_size = 100000

    sp = Species(species_id)
    midas_db = MIDAS_DB("design_chunks", "gtdb")

    tsprint("1")
    sp.get_repgenome(midas_db)
    tsprint("1")
    sp.compute_snps_chunks_run(midas_db, chunk_size)
    chunks_of_sites = load_chunks_cache(sp.chunks_of_sites_fp)
    tsprint("2")
    print(sp.num_of_snps_chunks)
    tsprint("2")
    print(chunks_of_sites[-1])

    tsprint("3")
    sp.get_repgenome(midas_db)
    tsprint("3")
    sp.compute_snps_chunks_merge(midas_db, chunk_size)
    tsprint("4")
    chunks_of_sites = load_chunks_cache(sp.chunks_of_sites_fp)
    tsprint("4")
    print(sp.num_of_snps_chunks)
    print(chunks_of_sites[-1])

    tsprint("5")
    sp.compute_gene_chunks(midas_db, chunk_size)
    tsprint("5")
    chunks_of_centroids = load_chunks_cache(sp.chunks_of_centroids_fp)
    tsprint("6")
    print(chunks_of_centroids[-1])
