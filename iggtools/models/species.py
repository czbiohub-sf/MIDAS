#!/usr/bin/env python3
from collections import defaultdict
from math import floor
import Bio.SeqIO

from iggtools.common.utils import tsprint, num_physical_cores, InputStream, OutputStream, command, select_from_tsv
from iggtools.models.uhgg import MIDAS_IGGDB
from iggtools.models.sample import Sample
from iggtools.params.schemas import genes_feature_schema


class Species:
    """ Base class for species """
    def __init__(self, species_id):
        self.id = species_id
        # SNPs chunks design related
        self.contigs = defaultdict(dict)
        self.chunks_of_sites = defaultdict(list)
        self.num_of_chunks = 0
        self.priority_chunks = []
        # MERGE related: select_species
        self.list_of_samples = [] # samples associate with current species
        self.samples_count = 0
        self.list_of_samples_depth = [] # mean_coverage
        # GENE Features for representative genomes
        self.genes_feature = None # indexed by contig_id
        self.genes_boundary = None # indexed by contig_id
        self.genes_sequence = None # indexed by gene_id


    def design_snps_chunks(self, midas_iggdb, chunk_size):
        """ Given the Genome and chunk_size, the structure of the chunks are the same.
            Each chunk is indexed by (species_id, chunk_id) """
        self.contigs = get_representative_genome(self.id, midas_iggdb)
        species_id = self.id

        # Start with full chunks
        chunk_id = 0
        dict_of_packed_args = defaultdict(list)
        priority_chunks = []
        unassigned_contigs = defaultdict(dict)
        for contig in self.contigs.values():
            contig_id = contig["contig_id"]
            contig_length = contig["contig_len"]
            #contig_seq = contig["contig_seq"]
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
                        dict_of_packed_args[chunk_id] = [(species_id, chunk_id, contig_id, ci, ci+chunk_size, count_flag)] #contig_seq[ci:ci+chunk_size]
                        if count_flag:
                            priority_chunks.append(chunk_id)
                        chunk_id += 1
        # Partition unassigned short contigs into subsets
        dict_of_chunks, chunk_id = partition_contigs_into_chunks(unassigned_contigs, chunk_size, chunk_id)

        # Add the partitioned subsets to chunks
        for chunk_dict in dict_of_chunks.values():
            _chunk_id = chunk_dict["chunk_id"]
            list_of_contigs = chunk_dict["contigs_id"]
            for _contig_id in list_of_contigs:
                cstart = unassigned_contigs[_contig_id]["contig_start"]
                cend = unassigned_contigs[_contig_id]["contig_end"]
                #cseq = self.contigs[_contig_id]["contig_seq"][cstart:cend]
                cflag = unassigned_contigs[_contig_id]["compute_reads"]
                dict_of_packed_args[_chunk_id].append((species_id, _chunk_id, _contig_id, cstart, cend, cflag)) #cseq
        assert chunk_id == _chunk_id+1

        # Finally the merge jobs
        dict_of_packed_args[-1] = (species_id, -1)

        self.num_of_chunks = chunk_id
        self.chunks_of_sites = dict_of_packed_args
        self.priority_chunks = priority_chunks

        return True


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
        self.genes_sequence = scan_genes(genes_seq_file)


    def prepare_annotation(self, genes_feature_file, genes_seq_file):
        self.fetch_genes_sequence(genes_seq_file)
        self.fetch_genes_feature(genes_feature_file)
        self.compute_genes_boundary()


def scan_contigs(contig_fasta_file):
    """ Scan the genome file to get contig_id and contig_seq as ref_seq """
    contigs = {}
    with InputStream(contig_fasta_file) as file:
        for rec in Bio.SeqIO.parse(file, 'fasta'):
            contigs[rec.id] = {
                "contig_id": rec.id,
                "contig_len": len(rec.seq),
                "contig_seq": str(rec.seq),
            }
    return contigs


def get_representative_genome(species_id, midas_iggdb):
    representative_genomes = midas_iggdb.fetch_files("prokka_genome", [species_id])
    contigs = scan_contigs(representative_genomes[species_id])
    return contigs


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


def scan_genes(fasta_file):
    """ Scan the genome file to get contig_id and contig_seq as ref_seq """
    genes = {}
    with InputStream(fasta_file) as file:
        for rec in Bio.SeqIO.parse(file, 'fasta'):
            genes[rec.id] = {
                "gene_id": rec.id,
                "gene_len": len(rec.seq),
                "gene_seq": str(rec.seq),
            }
    return genes


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
        gene_boundaries[contig_id] = {"genes": list(feature_ranges_sorted.keys()), "boundaries": boundaries} # "offsets": gene_offsets, "ranges": feature_ranges_sorted}
    return gene_boundaries
