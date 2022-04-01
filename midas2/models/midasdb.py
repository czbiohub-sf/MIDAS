#!/usr/bin/env python3
import os
from midas2.common.utils import download_reference, command, multithreading_map
from midas2.models.uhgg import UHGG
from midas2.params.inputs import MIDASDB_DICT, MARKER_FILE_EXTS, marker_set


# The "output" or built MIDAS DB layout in S3.
# See https://github.com/czbiohub/MIDAS2.0/wiki/MIDAS-DB#target-layout-in-s3
def get_midasdb_layout(species_id="", genome_id="", component=""):
    return {
        # Input: table of content and collections of genomes
        "table_of_contents":             f"genomes.tsv",
        "imported_genome":               f"cleaned_imports/{species_id}/{genome_id}/{genome_id}.{component}",

        "representative_genome":         f"gene_annotations/{species_id}/{genome_id}/{genome_id}.fna",
        "annotation_file":               f"gene_annotations/{species_id}/{genome_id}/{genome_id}.{component}",
        "annotation_fna":                f"gene_annotations/{species_id}/{genome_id}/{genome_id}.fna",
        "annotation_ffn":                f"gene_annotations/{species_id}/{genome_id}/{genome_id}.ffn",
        "annotation_genes":              f"gene_annotations/{species_id}/{genome_id}/{genome_id}.genes",
        "annotation_faa":                f"gene_annotations/{species_id}/{genome_id}/{genome_id}.faa",
        "annotation_gff":                f"gene_annotations/{species_id}/{genome_id}/{genome_id}.gff",
        "annotation_tsv":                f"gene_annotations/{species_id}/{genome_id}/{genome_id}.tsv",
        "annotation_log":                f"gene_annotations/{species_id}/{genome_id}/annotate_genome.log",
        "gene_features_log":             f"gene_annotations/{species_id}/{genome_id}/build_gene_features.log",

        "marker_db":                     [f"markers/{marker_set}/{marker_set}.{ext}" for ext in MARKER_FILE_EXTS],
        "marker_db_hmm":                 f"markers_models/{marker_set}/marker_genes.hmm",
        "marker_db_hmm_cutoffs":         f"markers_models/{marker_set}/marker_genes.mapping_cutoffs",
        "build_markerdb_log":            f"markers/{marker_set}/build_markerdb.log",

        "pangenome_file":                f"pangenomes/{species_id}/{component}",
        "pangenome_centroids":           f"pangenomes/{species_id}/centroids.ffn",
        "pangenome_cluster_info":        f"pangenomes/{species_id}/cluster_info.txt",
        "pangenome_log":                 f"pangenomes/{species_id}/pangenome_build.log",
        "pangenome_genes_info":          f"pangenomes/{species_id}/gene_info.txt",
        "pangenome_genes_len":           f"pangenomes/{species_id}/genes.len",
        "cluster_info_log":              f"pangenomes/{species_id}/cluster_info.log",

        "chunks_centroids":              f"chunks/genes/chunksize.{component}/{species_id}.json",
        "chunks_sites_run":              f"chunks/sites/run/chunksize.{component}/{species_id}/{genome_id}.json",
        "chunks_sites_merge":            f"chunks/sites/merge/chunksize.{component}/{species_id}/{genome_id}.json",

        "chunks_merge_list_of_contigs":  f"temp/chunks_merge_sites/chunksize.{component}/{species_id}/cid.{genome_id}_list_of_contigs",

        "marker_genes":                  f"markers/{marker_set}/temp/{species_id}/{genome_id}/{genome_id}.{component}",
        "marker_genes_seq":              f"markers/{marker_set}/temp/{species_id}/{genome_id}/{genome_id}.markers.fa",
        "marker_genes_map":              f"markers/{marker_set}/temp/{species_id}/{genome_id}/{genome_id}.markers.map",
        "marker_genes_hmmsearch":        f"markers/{marker_set}/temp/{species_id}/{genome_id}/{genome_id}.hmmsearch",
        "marker_genes_log":              f"markers/{marker_set}/temp/{species_id}/{genome_id}/build_marker_genes.log",
    }


class MIDAS_DB: # pylint: disable=too-few-public-methods

    def __init__(self, db_dir=".", db_name="uhgg", num_cores=1):
        self.db_dir = os.path.abspath(db_dir) # local dir
        self.db_name = MIDASDB_DICT[db_name] # s3 dir
        self.num_cores = num_cores
        self.local_toc = self.fetch_files("table_of_contents")

        self.uhgg = UHGG(self.local_toc)


    def get_repgenome_id(self, species_id):
        return self.uhgg.fetch_repgenome_id(species_id)
        #return self.uhgg.representatives[species_id]


    def construct_local_path(self, filename, species_id="", genome_id="", component=""):
        local_path = _get_local_path(get_midasdb_layout(species_id, genome_id, component)[filename], self.db_dir)
        return local_path


    def construct_dest_path(self, filename, species_id="", genome_id="", component=""):
        s3_path = _get_dest_path(get_midasdb_layout(species_id, genome_id, component)[filename], self.db_name)
        return s3_path


    def construct_file_tuple(self, filename, species_id="", genome_id="", component=""):
        local_path = self.construct_local_path(filename, species_id, genome_id, component)
        s3_path = self.construct_dest_path(filename, species_id, genome_id, component)
        return (s3_path, local_path)


    def fetch_file(self, filename, species_id="", genome_id="", component=""):
        (s3_path, local_path) = self.construct_file_tuple(filename, species_id, genome_id, component)
        return _fetch_file_from_s3((s3_path, local_path))


    def fetch_files(self, filename, list_of_species=""):
        """ Fetch files from S3 to local MIDAS_DB """
        if list_of_species:
            args_list = []
            for species_id in list_of_species:
                genome_id = self.get_repgenome_id(species_id)
                args_list.append(self.construct_file_tuple(filename, species_id, genome_id))

            if len(list_of_species) > 1:
                _fetched_files = multithreading_map(_fetch_file_from_s3, args_list, self.num_cores)
            else:
                _fetched_files = [_fetch_file_from_s3(args_list[0])]

            fetched_files = {}
            for species_index, species_id in enumerate(list_of_species):
                fetched_files[species_id] = _fetched_files[species_index]
            return fetched_files

        if filename == "marker_db":
            l1 = self.construct_dest_path(filename)
            l2 = self.construct_local_path(filename)
            list_of_file_tuples = list(zip(l1, l2))
            fetched_files = [_fetch_file_from_s3(_) for _ in list_of_file_tuples]
            return dict(zip(MARKER_FILE_EXTS, fetched_files))

        # Fetch single file with filename as the MIDAS DB layout key
        (s3_path, local_path) = self.construct_file_tuple(filename)
        return _fetch_file_from_s3((s3_path, local_path))


    def fetch_chunks(self, filename, list_of_species="", chunk_size=1000000):
        assert len(list_of_species) > 0, f"Empty list of species to download"
        args_list = []
        for species_id in list_of_species:
            genome_id = self.get_repgenome_id(species_id)
            args_list.append(self.construct_file_tuple(filename, species_id, genome_id, chunk_size))
        _fetched_files = multithreading_map(_fetch_file_from_s3, args_list, self.num_cores)

        fetched_files = {}
        for species_index, species_id in enumerate(list_of_species):
            fetched_files[species_id] = _fetched_files[species_index]
        return fetched_files



    def get_target_layout(self, filename, remote, species_id="", genome_id="", component=""):
        local_path = self.construct_local_path(filename, species_id, genome_id, component)
        s3_path = self.construct_dest_path(filename, species_id, genome_id, component)
        target_path = s3_path if remote else local_path
        return target_path


def _get_dest_path(file_name, db_name):
    """ Append S3 MIDAS_DB path to file(s) """
    if isinstance(file_name, list):
        return [os.path.join(db_name, f"{fn}.lz4") for fn in file_name]
    return os.path.join(db_name, f"{file_name}.lz4")


def _get_local_path(file_name, db_dir):
    """ Append local MIDAS_DB path to file(s) """
    if isinstance(file_name, list):
        return [os.path.join(db_dir, f"{fn}") for fn in file_name]
    return os.path.join(db_dir, f"{file_name}")


def _fetch_file_from_s3(packed_args):
    """ Fetch local path from AWS S3 dest path """
    s3_path, local_path = packed_args
    local_dir = os.path.dirname(local_path)
    if not os.path.isdir(local_dir):
        command(f"mkdir -p {local_dir}")
    if os.path.exists(local_path):
        return local_path
    return download_reference(s3_path, local_dir)
