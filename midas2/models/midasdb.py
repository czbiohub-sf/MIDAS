#!/usr/bin/env python3
import os
import json
from hashlib import md5
from collections import defaultdict

from midas2.common.utils import select_from_tsv, sorted_dict, InputStream, download_reference, command, multithreading_map, download_tarball
from midas2.params.inputs import MIDASDB_DICT, MARKER_FILE_EXTS, MD5SUM_JSON, marker_set
from midas2.params.outputs import genomes as TABLE_OF_CONTENTS


DEFAULT_SITE_CHUNK_SIZE = 1000000
DEFAULT_GENE_CHUNK_SIZE = 50000
DEFAULT_CHUNKS = [DEFAULT_GENE_CHUNK_SIZE, DEFAULT_SITE_CHUNK_SIZE, DEFAULT_SITE_CHUNK_SIZE]


# The Target Layout of MIDAS Reference Database both on S3 and locally.
def get_midasdb_layout(species_id="", genome_id="", component=""):
    return {
        # Input: table of content and collections of genomes
        "table_of_contents":             "genomes.tsv",
        "raw_genome":                    f"mags/{species_id}/{genome_id}.{component}",
        "imported_genome":               f"cleaned_imports/{species_id}/{genome_id}/{genome_id}.{component}",
        "import_log":                    f"cleaned_imports/{species_id}/{genome_id}/import_genome.log",

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
        "pangenome_tempfile":            f"pangenomes/{species_id}/temp/{genome_id}/{component}",
        "pangenome_log":                 f"pangenomes/{species_id}/pangenome_build.log",
        "recluster_log":                 f"pangenomes/{species_id}/recluster_centroids.log",

        "pangenome_marker_map":          f"pangenomes/{species_id}/temp/markers.map",
        "pangenome_centroids":           f"pangenomes/{species_id}/centroids.ffn",
        "pangenome_genes_info":          f"pangenomes/{species_id}/genes_info.tsv",
        "cluster_xx_info":               f"pangenomes/{species_id}/augment/clusters_{component}_info.tsv",
        "pangenome_contigs_len":         f"pangenomes/{species_id}/contigs.len",
        "pangenome_genes_annot":         f"pangenomes/{species_id}/genes_annotated.tsv",
        "cluster_xx_annot":              f"pangenomes/{species_id}/annotation/clusters_{component}_annot.tsv",
        "pangenome_cluster_xx":          f"pangenomes/{species_id}/clusters_{component}_info.tsv",

        "augment_log":                   f"pangenomes/{species_id}/augment.log",
        "contig_length_log":             f"pangenomes/{species_id}/contig_length.log",
        "enhance_log":                   f"pangenomes/{species_id}/pangenome_enhance.log",
        "prune_log":                     f"pangenomes/{species_id}/pruned/pangenome_prune.log",

        "genomad_virus_genes":           f"pangenomes_annotation/01_mge/{species_id}/{genome_id}/genomad_output/{genome_id}_summary/{genome_id}_virus_genes.tsv",
        "genomad_virus_summary":         f"pangenomes_annotation/01_mge/{species_id}/{genome_id}/genomad_output/{genome_id}_summary/{genome_id}_virus_summary.tsv",
        "genomad_plasmid_genes":         f"pangenomes_annotation/01_mge/{species_id}/{genome_id}/genomad_output/{genome_id}_summary/{genome_id}_plasmid_genes.tsv",
        "genomad_plasmid_summary":       f"pangenomes_annotation/01_mge/{species_id}/{genome_id}/genomad_output/{genome_id}_summary/{genome_id}_plasmid_summary.tsv",
        "mefinder_results":              f"pangenomes_annotation/01_mge/{species_id}/{genome_id}/mefinder_output/mefinder.csv",
        "resfinder_results":             f"pangenomes_annotation/01_mge/{species_id}/{genome_id}/resfinder_output/ResFinder_results_tab.txt",
        "eggnog_results":                f"pangenomes_annotation/02_eggnog/{species_id}/{species_id}.emapper.annotations",
        "panannot_tempfile":             f"pangenomes_annotation/03_processed/{species_id}/{component}/{genome_id}",

        "panannot_genomad_virus":        f"pangenomes/{species_id}/annotation/genomad_virus.tsv",
        "panannot_genomad_plasmid":      f"pangenomes/{species_id}/annotation/genomad_plasmid.tsv",
        "panannot_mefinder":             f"pangenomes/{species_id}/annotation/mefinder.tsv",
        "panannot_resfinder":            f"pangenomes/{species_id}/annotation/resfinder.tsv",
        "panannot_eggnog":               f"pangenomes/{species_id}/annotation/eggnog.tsv",
        "pruned_centroids":              f"pangenomes/{species_id}/pruned/centroids_by.{genome_id}.{component}.ffn",
        "pruned_centroids_rs":           f"pangenomes/{species_id}/pruned/centroids_by.{genome_id}.{component}.rmsig.ffn",

        "marker_genes":                  f"markers/{marker_set}/temp/{species_id}/{genome_id}/{genome_id}.{component}",
        "marker_genes_seq":              f"markers/{marker_set}/temp/{species_id}/{genome_id}/{genome_id}.markers.fa",
        "marker_genes_map":              f"markers/{marker_set}/temp/{species_id}/{genome_id}/{genome_id}.markers.map",
        "marker_genes_hmmsearch":        f"markers/{marker_set}/temp/{species_id}/{genome_id}/{genome_id}.hmmsearch",
        "marker_genes_log":              f"markers/{marker_set}/temp/{species_id}/{genome_id}/build_marker_genes.log",

        "chunks_sites_run":              f"chunks/sites/run/chunksize.{component}/{species_id}/{genome_id}.json",
        "chunks_sites_merge":            f"chunks/sites/merge/chunksize.{component}/{species_id}/{genome_id}.json",
        "chunks_contig_lists":           f"temp/chunksize.{component}/{species_id}/cid.{genome_id}_list_of_contigs",
    }



def get_tarball_layout(species_id="", genome_id=""): # target folder
    return {
        "toc":                           "genomes.tsv",
        "metadata":                      "metadata.tsv",
        "md5sum":                        "md5sum.json",
        "repgenome":                     f"gene_annotations/{species_id}/{genome_id}",
        "pangenome":                     f"pangenomes/{species_id}",
        "markerdb":                      f"markers/{marker_set}",
        "markerdb_models":               f"markers_models/{marker_set}",
        "chunks":                        "chunks"
    }


# TODO: update the pangenome componenets
tarball_mapping = {
    "toc":                           ["table_of_contents"],
    "repgenome":                     ["annotation_fna", "annotation_ffn", "annotation_genes"],
    "pangenome":                     ["pangenome_centroids", "pangenome_genes_info"],
    "markerdb":                      MARKER_FILE_EXTS,
    "markerdb_models":               ["hmm", "hmm_cutoffs"],
    "chunks":                        ["chunks_centroids", "chunks_sites_run", "chunks_sites_merge"],
}


class UHGG:  # pylint: disable=too-few-public-methods
    """
    A model for the external prokaryotic collection of genomes.
    """
    def __init__(self, table_of_contents_tsv="", dbname="uhgg"):
        if table_of_contents_tsv:
            self.toc_tsv = table_of_contents_tsv # Read in local genomes.tsv
        else:
            self.toc_tsv = TABLE_OF_CONTENTS(dbname) # Fetch from S3 only

        self.species, self.representatives, self.genomes = _UHGG_load(self.toc_tsv)

    def fetch_repgenome_id(self, species_id):
        return self.representatives[species_id]


class MIDAS_DB: # pylint: disable=too-few-public-methods

    def __init__(self, db_dir=".", db_name="uhgg", num_cores=1):
        self.db_dir = os.path.abspath(db_dir) # local dir
        self.db_name = MIDASDB_DICT[db_name] # remote dir
        self.num_cores = num_cores
        self.has_md5sum = db_name in MD5SUM_JSON
        self.md5sum = load_json(self.fetch_files("md5sum")) if self.has_md5sum else None
        self.local_toc = self.fetch_files("toc" if self.has_md5sum else "table_of_contents")
        self.uhgg = UHGG(self.local_toc, db_name)


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
        """ Fetch Single File Only From *S3* Bucket. And return local single-file path. """
        if species_id and not genome_id:
            genome_id = self.get_repgenome_id(species_id)
        (s3_path, local_path) = self.construct_file_tuple(filename, species_id, genome_id, component)
        return _fetch_file_from_s3((s3_path, local_path))


    def fetch_files(self, filename, list_of_species="", rep_only=True):
        """ Two purposes: fetch remote file (if not exists) and return a dictionary. """
        if filename in get_tarball_layout():
            # Return Fetched Directory
            return self.fetch_tarball(filename, list_of_species)
        # Return Fetched File(s)
        return self.fetch_individual_files(filename, list_of_species, rep_only)


    def fetch_tarball(self, filename, list_of_species):
        """ Fetch tarball from GIDB server and Check MD5SUM """
        # Rep-genome and Pan-genome
        if list_of_species:
            args_list = []
            for species_id in list_of_species:
                genome_id = self.get_repgenome_id(species_id)
                args_list.append(self.construct_file_tuple(filename, species_id, genome_id))
            if len(args_list) > 1:
                _fetched_files = multithreading_map(_fetch_file_from_s3, args_list, self.num_cores)
            else:
                _fetched_files = [_fetch_file_from_s3(args_list[0])]

            fetched_files = dict(zip(list_of_species, _fetched_files))
            if self.has_md5sum:
                for species_id in list_of_species:
                    fetched_filenames = tarball_mapping[filename]
                    for _filename in fetched_filenames:
                        genome_id = self.get_repgenome_id(species_id)
                        _fetched_file = self.get_target_layout(_filename, False, species_id, genome_id)
                        md5_fetched = file_md5sum(_fetched_file)
                        md5_lookup = self.md5sum[filename][species_id][_filename]
                        assert md5_fetched == md5_lookup, f"Error for downloading {_fetched_file} from {filename}. Please delete the folder and redownload."
            return fetched_files

        # Single Copy Marker Genes DB
        if filename == "markerdb":
            fetched_dir = _fetch_file_from_s3(self.construct_file_tuple(filename))
            fetched_files = self.get_target_layout("marker_db", False)
            fetched_files = dict(zip(MARKER_FILE_EXTS, fetched_files))
            if self.has_md5sum:
                for _ in MARKER_FILE_EXTS:
                    _fetched_file = fetched_files[_]
                    md5_fetched = file_md5sum(_fetched_file)
                    md5_lookup = self.md5sum[filename][_]
                    assert md5_fetched == md5_lookup, f"Error for downloadding {_fetched_file} from {filename}. Please delete the folder and redownload."
            return {filename: fetched_dir}

        if filename == "markerdb_models":
            fetched_dir = _fetch_file_from_s3(self.construct_file_tuple(filename))
            if self.has_md5sum:
                for _ in tarball_mapping[filename]:
                    _fetched_file = self.get_target_layout(f"marker_db_{_}", False)
                    md5_fetched = file_md5sum(_fetched_file)
                    md5_lookup = self.md5sum[filename][_]
                    assert md5_fetched == md5_lookup, f"Error for downloadding {_fetched_file} from {filename}. Please delete the folder and redownload."
            return {filename: fetched_dir}

        # Chunks
        if filename == "chunks":
            fetched_dir = _fetch_file_from_s3(self.construct_file_tuple(filename))
            if self.has_md5sum:
                fetched_filenames = tarball_mapping[filename]
                rep_genomes = self.uhgg.representatives
                for sid, gid in rep_genomes.items():
                    for i, ct in enumerate(fetched_filenames):
                        _fetched_file = self.get_target_layout(ct, False, sid, gid, DEFAULT_CHUNKS[i])
                        md5_fetched = file_md5sum(_fetched_file)
                        md5_lookup = self.md5sum[filename][sid][ct]
                        assert md5_fetched == md5_lookup, f"Error for downloadding {_fetched_file} from {filename}. Please delete the folder and redownload."
            return {filename: fetched_dir}

        # Single File: key of the tarball layout
        _fetched_file = _fetch_file_from_s3(self.construct_file_tuple(filename))
        md5_fetched = file_md5sum(_fetched_file)
        if filename != "md5sum":
            md5_lookup = self.md5sum[filename]
        else:
            md5_lookup = MD5SUM_JSON[os.path.basename(self.db_name)]
        assert md5_fetched == md5_lookup, f"Error for downloadding {_fetched_file} from {filename}. Please delete the folder and redownload."
        return _fetched_file


    def fetch_individual_files(self, filename, list_of_species, rep_only):
        """ Fetch S3 individual file (if not exists) and/or return local file path in a dictionary. """
        if list_of_species:
            if not rep_only:
                assert len(list_of_species) == 1, "Only download all genomes for single species"

            args_list = []
            for species_id in list_of_species:
                if rep_only:
                    genome_id = self.get_repgenome_id(species_id)
                    args_list.append(self.construct_file_tuple(filename, species_id, genome_id))
                else:
                    # During pan-genome and marker-db built, collect all the genomes for given species.
                    # Either from S3 or locally, therefore no need to check md5sum.
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
    if db_name.startswith("https://"):
        compress_cmd = "tar.gz"
    else:
        compress_cmd = "lz4"
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
    # local_path: single-file or single-folder
    # local_dir: parent directory of local_path
    # s3_path: https or s3 link of tarball of LZ4 file
    if not os.path.isdir(local_dir):
        command(f"mkdir -p {local_dir}")
    if os.path.exists(local_path):
        return local_path
    if s3_path.startswith("s3://"):
        return download_reference(s3_path, local_dir)
    return download_tarball(s3_path, local_dir)


def load_json(jsonfile):
    with InputStream(jsonfile) as stream:
        md5sum_dict = json.load(stream)
    return md5sum_dict


def file_md5sum(local_file):
    return md5(open(local_file, "rb").read()).hexdigest()
