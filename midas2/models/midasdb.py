#!/usr/bin/env python3
import os
from collections import defaultdict
from midas2.common.utils import select_from_tsv, sorted_dict, InputStream, download_reference, command, multithreading_map, download_tarball
from midas2.params.inputs import MIDASDB_DICT, MARKER_FILE_EXTS, marker_set
from midas2.params.outputs import genomes as TABLE_OF_CONTENTS


# The Target Layout of MIDAS Reference Database both on S3 and locally.
def get_midasdb_layout(species_id="", genome_id="", component=""):
    return {
        # Input: table of content and collections of genomes
        "table_of_contents":             f"genomes.tsv",
        "imported_genome":               f"cleaned_imports/{species_id}/{genome_id}/{genome_id}.{component}",

        "annotation_file":               f"gene_annotations/{species_id}/{genome_id}/{genome_id}.{component}",
        "representative_genome":         f"gene_annotations/{species_id}/{genome_id}/{genome_id}.fna",
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

        "marker_genes":                  f"markers/{marker_set}/temp/{species_id}/{genome_id}/{genome_id}.{component}",
        "marker_genes_seq":              f"markers/{marker_set}/temp/{species_id}/{genome_id}/{genome_id}.markers.fa",
        "marker_genes_map":              f"markers/{marker_set}/temp/{species_id}/{genome_id}/{genome_id}.markers.map",
        "marker_genes_hmmsearch":        f"markers/{marker_set}/temp/{species_id}/{genome_id}/{genome_id}.hmmsearch",
        "marker_genes_log":              f"markers/{marker_set}/temp/{species_id}/{genome_id}/build_marker_genes.log",

        "chunks_centroids":              f"chunks/genes/chunksize.{component}/{species_id}.json",
        "chunks_sites_run":              f"chunks/sites/run/chunksize.{component}/{species_id}/{genome_id}.json",
        "chunks_sites_merge":            f"chunks/sites/merge/chunksize.{component}/{species_id}/{genome_id}.json",
        "chunks_contig_lists":           f"temp/chunksize.{component}/{species_id}/cid.{genome_id}_list_of_contigs",
    }


def get_tarball_layout(species_id="", genome_id=""):
    return {
        # Target Folder
        "table_of_contents":             f"genomes.tsv",
        "repgenome":                     f"gene_annotations/{species_id}/{genome_id}",
        "pangenome":                     f"pangenomes/{species_id}",
        "markerdb":                      f"markers/{marker_set}",
        "markerdb_models":               f"markers_models/{marker_set}",
        "chunks":                        f"chunks"
    }


class UHGG:  # pylint: disable=too-few-public-methods
    """
    A model for the external prokaryotic collection of genomes.
    """
    def __init__(self, table_of_contents_tsv=""):
        if table_of_contents_tsv:
            self.toc_tsv = table_of_contents_tsv # Read in local genomes.tsv
        else:
            self.toc_tsv = TABLE_OF_CONTENTS # Fetch from S3

        self.species, self.representatives, self.genomes = _UHGG_load(self.toc_tsv)

    def fetch_repgenome_id(self, species_id):
        return self.representatives[species_id]


class MIDAS_DB: # pylint: disable=too-few-public-methods

    def __init__(self, db_dir=".", db_name="uhgg", num_cores=1):
        self.db_dir = os.path.abspath(db_dir) # local dir
        self.db_name = MIDASDB_DICT[db_name] # remote dir
        self.num_cores = num_cores
        self.local_toc = self.fetch_files("table_of_contents")
        self.uhgg = UHGG(self.local_toc)


    def get_repgenome_id(self, species_id):
        return self.uhgg.fetch_repgenome_id(species_id)


    def construct_local_path(self, filename, species_id="", genome_id="", component=""):
        # We need single-file access for server download files as well
        if filename in get_tarball_layout():
            local_path = _get_local_path(get_tarball_layout(species_id, genome_id)[filename], self.db_dir)
        else:
            local_path = _get_local_path(get_midasdb_layout(species_id, genome_id, component)[filename], self.db_dir)
        return local_path


    def construct_dest_path(self, filename, species_id="", genome_id="", component=""):
        if filename in get_tarball_layout():
            remote_path = _get_dest_path(get_tarball_layout(species_id, genome_id)[filename], self.db_name)
        else:
            remote_path = _get_dest_path(get_midasdb_layout(species_id, genome_id, component)[filename], self.db_name)
        return remote_path


    def get_target_layout(self, filename, remote, species_id="", genome_id="", component=""):
        local_path = self.construct_local_path(filename, species_id, genome_id, component)
        if not remote:
            return local_path

        s3_path = self.construct_dest_path(filename, species_id, genome_id, component)
        return s3_path


    def construct_file_tuple(self, filename, species_id="", genome_id="", component=""):
        local_path = self.construct_local_path(filename, species_id, genome_id, component)
        s3_path = self.construct_dest_path(filename, species_id, genome_id, component)
        return (s3_path, local_path)


    def fetch_file(self, filename, species_id="", genome_id="", component=""):
        """ Fetch single file only from S3 """
        (s3_path, local_path) = self.construct_file_tuple(filename, species_id, genome_id, component)
        return _fetch_file_from_s3((s3_path, local_path))


    def fetch_files(self, filename, list_of_species="", rep_only=True):
        """ Fetch files for list of species """
        if list_of_species:
            if not rep_only:
                assert len(list_of_species) == 1, f"Only download all genomes for single species"

            args_list = []
            for species_id in list_of_species:
                if rep_only:
                    genome_id = self.get_repgenome_id(species_id)
                    args_list.append(self.construct_file_tuple(filename, species_id, genome_id))
                else:
                    species_genomes_ids = self.uhgg.species[species_id].keys()
                    args_list.extend([self.construct_file_tuple(filename, species_id, genome_id) for genome_id in species_genomes_ids])

            if len(args_list) > 1:
                _fetched_files = multithreading_map(_fetch_file_from_s3, args_list, self.num_cores)
            else:
                _fetched_files = [_fetch_file_from_s3(args_list[0])]

            if rep_only:
                fetched_files = dict(zip(list_of_species, _fetched_files))
            else:
                fetched_files = dict(zip(species_genomes_ids, _fetched_files))
            return fetched_files

        if filename == "marker_db":
            l1 = self.construct_dest_path(filename)
            l2 = self.construct_local_path(filename)
            list_of_file_tuples = list(zip(l1, l2))
            fetched_files = [_fetch_file_from_s3(_) for _ in list_of_file_tuples]
            return dict(zip(MARKER_FILE_EXTS, fetched_files))

        # Fetch single file with filename as the MIDASDB Target Layout key
        return _fetch_file_from_s3(self.construct_file_tuple(filename))


def _UHGG_load(toc_tsv, deep_sort=False):
    species = defaultdict(dict)
    representatives = {}
    genomes = {}
    with InputStream(toc_tsv) as table_of_contents:
        for row in select_from_tsv(table_of_contents, selected_columns=["genome", "species", "representative", "genome_is_representative"]):
            genome_id, species_id, representative_id, _ = row
            species[species_id][genome_id] = row
            representatives[species_id] = representative_id
            genomes[genome_id] = species_id
    if deep_sort:
        for sid in species.keys():
            species[sid] = sorted_dict(species[sid])
        species = sorted_dict(species)
    return species, representatives, genomes


def _get_dest_path(file_name, db_name):
    """ Append S3 MIDAS_DB path to file(s) """
    if db_name.startswith("s3://"):
        compress_cmd = "lz4"
    else:
        compress_cmd = "tar.gz"

    if isinstance(file_name, list):
        return [os.path.join(db_name, f"{fn}.{compress_cmd}") for fn in file_name]
    return os.path.join(db_name, f"{file_name}.{compress_cmd}")


def _get_local_path(file_name, db_dir):
    """ Append local MIDAS_DB path to file(s) """
    if isinstance(file_name, list):
        return [os.path.join(db_dir, f"{fn}") for fn in file_name]
    return os.path.join(db_dir, f"{file_name}")


def _fetch_file_from_s3(packed_args):
    """ Fetch local path from Remote. Return local_path if already exist. """
    s3_path, local_path = packed_args
    local_dir = os.path.dirname(local_path)

    print(f"_fetch_file_from_s3: {s3_path}")
    print(f"_fetch_file_from_s3: {local_dir}")
    print(f"_fetch_file_from_s3: {local_path}")

    if not os.path.isdir(local_dir):
        command(f"mkdir -p {local_dir}")
    if os.path.exists(local_path):
        return local_path

    if s3_path.startswith("s3://"):
        return download_reference(s3_path, local_dir)
    return download_tarball(s3_path, local_dir)
