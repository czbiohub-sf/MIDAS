#!/usr/bin/env python3
import os
import json
from collections import defaultdict
from math import floor
import Bio.SeqIO

from iggtools.common.utils import tsprint, InputStream, OutputStream, command, select_from_tsv, cat_files
from iggtools.models.sample import Sample
from iggtools.params.schemas import genes_feature_schema, CLUSTER_INFO_SCHEMA


class Species:
    """ Base class for species """
    def __init__(self, species_id, cache=False):
        self.id = species_id
        self.cache = cache
        # SNPs chunks design related
        self.contigs = defaultdict(dict)
        self.max_contig_length = 0
        self.chunks_of_sites = defaultdict(list)
        self.num_of_sites_chunks = 0
        self.priority_chunks = []
        # MERGE related: select_species
        self.list_of_samples = [] # samples associate with current species
        self.samples_count = 0
        self.list_of_samples_depth = [] # mean_coverage
        # GENE Features for representative genomes
        self.genes_feature = None # indexed by contig_id
        self.genes_boundary = None # indexed by contig_id
        self.genes_sequence = None # indexed by gene_id
        # Pan genes
        #self.centroids = defaultdict(dict)
        self.cluster_info = None
        self.chunks_of_centroids = defaultdict(dict)
        self.num_of_genes_chunks = 0
        self.dict_of_genes_are_markers = []
        self.list_of_markers = []


    def design_genes_chunks(self, midas_db, chunk_size):
        """ Each chunk is indexed by (species_id, chunk_id) """
        species_id = self.id

        chunk_cache_fp = midas_db.get_target_layout("cache_gene_chunks", False, "chunks_of_centroids", species_id, chunk_size)
        if self.cache and os.path.exists(chunk_cache_fp):
            with InputStream(chunk_cache_fp) as stream:
                chunks_of_centroids = json.load(stream)
                chunk_id = len(chunks_of_centroids.keys())
        else:
            self.get_cluster_info(midas_db)
            genes_counter = 0
            chunk_id = 0

            curr_chunk_of_genes = dict()
            chunks_of_centroids = defaultdict(dict)
            for row in self.cluster_info.values():
                if not chunk_id*chunk_size <= genes_counter < (chunk_id+1)*chunk_size:
                    chunks_of_centroids[chunk_id] = curr_chunk_of_genes
                    chunk_id += 1
                    curr_chunk_of_genes = defaultdict()

                curr_chunk_of_genes[row["centroid_99"]] = row["centroid_99_length"]
                genes_counter += 1

            # Last chunk of centroids
            chunks_of_centroids[chunk_id] = curr_chunk_of_genes
            chunk_id += 1

            if self.cache:
                with OutputStream(chunk_cache_fp) as stream:
                    json.dump(chunks_of_centroids, stream)

        self.num_of_genes_chunks = chunk_id
        self.chunks_of_centroids = chunks_of_centroids

        return True


    def design_snps_chunks(self, midas_db, chunk_size):
        """ Given the Genome and chunk_size, the structure of the chunks are the same.
            Each chunk is indexed by (species_id, chunk_id) """

        species_id = self.id
        self.get_representative_genome(midas_db)

        # Start with full chunks
        chunk_id = 0
        dict_of_packed_args = defaultdict(list)
        priority_chunks = []
        unassigned_contigs = defaultdict(dict)
        max_contig_length = 0

        for contig in self.contigs.values():
            contig_id = contig["id"]
            contig_length = contig["length"]

            if contig_length > max_contig_length:
                max_contig_length = contig_length

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
                        dict_of_packed_args[chunk_id] = [(species_id, chunk_id, contig_id, ci, ci+chunk_size, count_flag, 0)]
                        if count_flag:
                            priority_chunks.append(chunk_id)
                        chunk_id += 1

        # Partition unassigned short contigs into subsets
        dict_of_chunks, chunk_id = partition_contigs_into_chunks(unassigned_contigs, chunk_size, chunk_id)

        # Add the partitioned subsets to chunks
        for chunk_dict in dict_of_chunks.values():
            _chunk_id = chunk_dict["chunk_id"]
            list_of_contigs = chunk_dict["contigs_id"]
            for wc_cidx, _cid in enumerate(list_of_contigs):
                cstart = unassigned_contigs[_cid]["contig_start"]
                cend = unassigned_contigs[_cid]["contig_end"]
                cflag = unassigned_contigs[_cid]["compute_reads"]
                dict_of_packed_args[_chunk_id].append((species_id, _chunk_id, _cid, cstart, cend, cflag, wc_cidx))
        assert chunk_id == _chunk_id+1

        # Finally the merge jobs
        dict_of_packed_args[-1] = (species_id, -1)

        self.num_of_sites_chunks = chunk_id
        self.chunks_of_sites = dict_of_packed_args
        self.priority_chunks = priority_chunks
        self.max_contig_length = max_contig_length

        return True


    def get_representative_genome(self, midas_db):
        species_id = self.id
        self.contigs = scan_fasta(midas_db.fetch_files("prokka_genome", [species_id])[species_id])


    def get_centroids(self, midas_db):
        species_id = self.id
        self.centroids = scan_fasta(midas_db.fetch_files("centroids", [species_id])[species_id])


    def get_cluster_info(self, midas_db):
        species_id = self.id
        cluster_info_path = midas_db.fetch_files("cluster_info", [species_id])[species_id]
        cluster_info = dict()
        with InputStream(cluster_info_path) as stream:
            for r in select_from_tsv(stream, selected_columns=CLUSTER_INFO_SCHEMA, result_structure=dict):
                cluster_info[r["centroid_99"]] = r
        self.cluster_info = cluster_info


    def fetch_samples_names(self):
        return [sample.sample_name for sample in self.list_of_samples]


    def fetch_samples_depth(self):
        self.list_of_samples_depth = [sample.profile[self.id]["mean_coverage"] for sample in self.list_of_samples]


    def fetch_genes_feature(self, features_file):
        self.genes_feature = read_gene_features(features_file)


    def compute_genes_boundary(self):
        assert self.genes_feature is not None, f"Need to fetch genes feature first"
        self.genes_boundary = generate_boundaries(self.genes_feature)


    def fetch_genes_sequence(self, genes_seq_file):
        self.genes_sequence = scan_fasta(genes_seq_file)


    def prepare_annotation(self, genes_feature_file, genes_seq_file):
        self.fetch_genes_sequence(genes_seq_file)
        self.fetch_genes_feature(genes_feature_file)
        self.compute_genes_boundary()


    def fetch_genes_are_markers(self, midas_db, chunk_size):
        species_id = self.id
        marker_centroids_file = midas_db.fetch_files("marker_centroids", [species_id])[species_id]

        chunk_cache_fp = midas_db.get_target_layout("cache_gene_chunks", False, "dict_of_genes_are_markers", species_id, chunk_size)
        if self.cache and os.path.exists(chunk_cache_fp):
            with InputStream(chunk_cache_fp) as stream:
                dict_of_genes_are_markers = json.load(stream)
                list_of_marker_genes = []
                for _ in dict_of_genes_are_markers.values():
                    if _["marker_id"] not in list_of_marker_genes:
                        list_of_marker_genes.append(_["marker_id"])
        else:
            dict_of_genes_are_markers = defaultdict(dict)
            list_of_marker_genes = []
            with InputStream(marker_centroids_file) as stream:
                for r in select_from_tsv(stream, selected_columns=["centroid_99", "marker_id"], result_structure=dict):
                    dict_of_genes_are_markers[r["centroid_99"]] = r
                    if r["marker_id"] not in list_of_marker_genes:
                        list_of_marker_genes.append(r["marker_id"])
            if self.cache:
                with OutputStream(chunk_cache_fp) as stream:
                    json.dump(dict_of_genes_are_markers, stream)

        self.dict_of_genes_are_markers = dict_of_genes_are_markers
        self.list_of_markers = list_of_marker_genes


def scan_fasta(fasta_file):
    """ Scan FASTA FILE to get seq and len """
    seqs = {}
    with InputStream(fasta_file) as file:
        for rec in Bio.SeqIO.parse(file, 'fasta'):
            seqs[rec.id] = {
                "id": rec.id,
                "length": len(rec.seq),
                "seq": str(rec.seq),
            }
    return seqs


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


def sort_list_of_species(list_of_species):
    """ Sort list_of_species by samples_count in descending order """
    species_sorted = sorted(((sp, len(sp.list_of_samples)) for sp in list_of_species), key=lambda x: x[1], reverse=True)
    return [sp[0] for sp in species_sorted]


def partition_contigs_into_chunks(unassigned_contigs, chunk_size, chunk_id):
    """ Similar to the problem of partition to K equal sum subsets """
    # Sort contigs by descending order of contig length
    sorted_contigs = {cid:cc["contig_length"] for cid, cc in sorted(unassigned_contigs.items(), key=lambda x: x[1]["contig_length"], reverse=True)}
    list_of_contigs_id = list(sorted_contigs.keys())
    list_of_contigs_length = list(sorted_contigs.values())
    # Initialize two pointers
    istart = 0
    jstart = len(list_of_contigs_length)-1
    prev_jstart = jstart + 1
    dict_of_chunks = defaultdict(dict)
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
        dict_of_chunks[chunk_id] = {"chunk_id": chunk_id, "contigs_id": curr_cids, "chunk_length": curr_chunk_length, "list_of_contigs_length": curr_clens}
        # Update the pointer and chunk_id
        list_of_partitoned_contigs = list_of_partitoned_contigs + curr_clens
        istart = istart + 1
        prev_jstart = jstart + 1
        chunk_id += 1
    assert len(list_of_partitoned_contigs) == len(list_of_contigs_length)
    assert set(list_of_partitoned_contigs) == set(list_of_contigs_length)
    return (dict_of_chunks, chunk_id)


def cluster_small_contigs(unassigned_contigs, chunk_size):
    """ Sort the contigs in ascending order and cluster small ones. The benefit is larger contigs mostly have individual chunks"""
    sorted_contigs = {cid:clen for cid, clen in sorted(unassigned_contigs.items(), key=lambda x: x[1])}
    list_of_contigs_id = list(sorted_contigs.keys())
    list_of_contigs_length = list(sorted_contigs.values())
    istart = 0
    prev_istart = 0
    chunk_id = 0
    dict_of_chunks = defaultdict(dict)
    t = []
    while istart < len(list_of_contigs_length):
        while sum(list_of_contigs_length[prev_istart:istart+1]) <= chunk_size and istart < len(list_of_contigs_length):
            istart += 1
        chunk_id += 1
        dict_of_chunks[chunk_id] = {"chunk_id": chunk_id, "contigs_id": list_of_contigs_id[prev_istart:istart],
                                    "chunk_length": sum(list_of_contigs_length[prev_istart:istart]),
                                    "contigs_len": list_of_contigs_length[prev_istart:istart]}
        t = t + list_of_contigs_length[prev_istart:istart]
        prev_istart = istart
    assert len(t) == len(list_of_contigs_length)
    assert set(t) == set(list_of_contigs_length)
    return dict_of_chunks


def read_gene_features(features_file):
    """ Read TAB-delimited *.genes files from gene_annotations """
    features = defaultdict(dict)
    with InputStream(features_file) as stream:
        for r in select_from_tsv(stream, selected_columns=genes_feature_schema, result_structure=dict):
            if r['gene_type'] == "CDS":
                features[r['contig_id']][r['gene_id']] = r ## gnl|Prokka|
    return features


def generate_boundaries(features):
    """ Given list of gene ranges, generate the desired, half-open boundaries by binary search """
    gene_boundaries = defaultdict(dict)
    for contig_id in features.keys():
        features_per_contig = features[contig_id]
        # Sort gene features by starting genomic position
        dict_of_feature_tuples = {fid: (feat['start'], feat['end']) for fid, feat in features_per_contig.items()}
        feature_ranges_sorted = dict(sorted(dict_of_feature_tuples.items(), key=lambda x: x[1][0], reverse=False))
        # For non-overlapping gene ranges, linear search would report the first met range (smallest)
        # therefore, we update the gene feature files for the genes with overlapping ranges with the adjacent genes before
        prev_gene = None
        gc = 0
        gene_offsets = dict()
        for gid, grange in feature_ranges_sorted.items():
            gc += 1
            if prev_gene is None:
                prev_gene = gid
                continue
            curr_gene = gid
            curr_start = feature_ranges_sorted[curr_gene][0]
            prev_end = feature_ranges_sorted[prev_gene][1]
            if curr_start <= prev_end:
                new_curr_start = prev_end + 1
                feature_ranges_sorted[curr_gene] = (new_curr_start, feature_ranges_sorted[curr_gene][1])
                gene_offsets[curr_gene] = new_curr_start - curr_start + 1
            prev_gene = curr_gene
        # Now, we are sure the sorted feature ranges don't have overlapping anymore
        feature_ranges_flat = tuple(_ for rt in tuple(feature_ranges_sorted.values()) for _ in rt)
        # Convert ranges into half-open intervals.
        boundaries = tuple(gr + 1 if idx%2 == 1 else gr for idx, gr in enumerate(feature_ranges_flat))
        gene_boundaries[contig_id] = {"genes": list(feature_ranges_sorted.keys()), "boundaries": boundaries}
    return gene_boundaries


def collect_units_per_chunk(sample, contig_counts_per_chunk, species_id, chunk_id, filename):
    ## Clean up sliced chunk temporary files
    headerless_sliced_path = sample.get_target_layout(filename, species_id, chunk_id)

    if contig_counts_per_chunk > 1:
        list_of_sliced_files = [sample.get_target_layout(f"{filename}_perc", species_id, chunk_id, cidx) for cidx in range(0, contig_counts_per_chunk)]
        cat_files(list_of_sliced_files, headerless_sliced_path, 20)
        for s_file in list_of_sliced_files:
            command(f"rm -rf {s_file}", quiet=True)
    else:
        sliced_file = sample.get_target_layout(f"{filename}_perc", species_id, chunk_id, 0)
        command(f"mv {sliced_file} {headerless_sliced_path}", quiet=True)
