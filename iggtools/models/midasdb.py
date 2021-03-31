#!/usr/bin/env python3
# A model for the UHGG collection of genomes (aka UHGG database).
import os
from iggtools.common.utils import download_reference, multithreading_map, command
from iggtools.params import outputs
from iggtools.models.uhgg import UHGG, get_uhgg_layout, destpath
from iggtools.params.schemas import MARKER_FILE_EXTS


class MIDAS_DB: # pylint: disable=too-few-public-methods

    def __init__(self, midas_db_dir=".", num_cores=1):
        self.midas_db_dir = midas_db_dir
        self.num_cores = num_cores
        self.local_toc = _fetch_file_from_s3((outputs.genomes, self.get_target_layout("genomes_toc", False)))
        self.uhgg = UHGG(self.local_toc)


    def get_target_layout(self, filename, remote, component="", species_id="", genome_id=""):
        local_path = get_uhgg_layout(species_id, component, genome_id)[filename]
        target_path = destpath(local_path) if remote else os.path.join(self.midas_db_dir, local_path)
        return target_path


    def fetch_files(self, filetype, list_of_species_ids=""):
        """ Fetch files from {inputs.igg} to local midas _db """
        args_list = []
        fetched_files = {}

        if filetype == "marker_db":
            for ext in MARKER_FILE_EXTS:
                s3_file = self.get_target_layout("marker_db", remote=True, component=ext)
                dest_file = self.get_target_layout("marker_db", remote=False, component=ext)
                fetched_files[ext] = _fetch_file_from_s3((s3_file, dest_file))
            return fetched_files

        # Download per species files
        if list_of_species_ids:
            for species_id in list_of_species_ids:
                if filetype == "contigs":
                    s3_file = self.get_target_layout("imported_genome_file", True, "fna", species_id, self.uhgg.representatives[species_id])
                    dest_file = self.get_target_layout("imported_genome_file", False, "fna", species_id, self.uhgg.representatives[species_id])
                if filetype == "centroids":
                    s3_file = self.get_target_layout("pangenome_file", True, "centroids.ffn", species_id)
                    dest_file = self.get_target_layout("pangenome_file", False, "centroids.ffn", species_id)
                if filetype == "genes_info":
                    s3_file = self.get_target_layout("pangenome_file", True, "gene_info.txt", species_id)
                    dest_file = self.get_target_layout("pangenome_file", False, "gene_info.txt", species_id)
                if filetype == "cluster_info":
                    s3_file = self.get_target_layout("pangenome_file", True, "cluster_info.txt", species_id)
                    dest_file = self.get_target_layout("pangenome_file", False, "cluster_info.txt", species_id)

                if filetype == "marker_genes_fa":
                    s3_file = self.get_target_layout("marker_genes", True, "markers.fa", species_id, self.uhgg.representatives[species_id])
                    dest_file = self.get_target_layout("marker_genes", False, "markers.fa", species_id, self.uhgg.representatives[species_id])
                if filetype == "marker_genes_map":
                    s3_file = self.get_target_layout("marker_genes", True, "markers.map", species_id, self.uhgg.representatives[species_id])
                    dest_file = self.get_target_layout("marker_genes", False, "markers.map", species_id, self.uhgg.representatives[species_id])

                if filetype == "marker_centroids":
                    s3_file = self.get_target_layout("marker_centroids", True, "txt", species_id)
                    dest_file = self.get_target_layout("marker_centroids", False, "txt", species_id)

                if filetype == "gene_feature":
                    s3_file = self.get_target_layout("annotation_file", True, "genes", species_id, self.uhgg.representatives[species_id])
                    dest_file = self.get_target_layout("annotation_file", False, "genes", species_id, self.uhgg.representatives[species_id])
                if filetype == "gene_seq":
                    s3_file = self.get_target_layout("annotation_file", True, "ffn", species_id, self.uhgg.representatives[species_id])
                    dest_file = self.get_target_layout("annotation_file", False, "ffn", species_id, self.uhgg.representatives[species_id])

                if filetype == "prokka_genome":
                    s3_file = self.get_target_layout("annotation_file", True, "fna", species_id, self.uhgg.representatives[species_id])
                    dest_file = self.get_target_layout("annotation_file", False, "fna", species_id, self.uhgg.representatives[species_id])

                args_list.append((s3_file, dest_file))

            _fetched_files = multithreading_map(_fetch_file_from_s3, args_list, num_threads=self.num_cores)
            for species_index, species_id in enumerate(list_of_species_ids):
                fetched_files[species_id] = _fetched_files[species_index]
            return fetched_files

        # Fetch single file with filetype as the uhgg layout key
        if filetype in get_uhgg_layout(species_id=""):
            s3_file = self.get_target_layout(filetype, True)
            dest_file = self.get_target_layout(filetype, False)
        else:
            s3_file = destpath(filetype)
            dest_file = os.path.join(self.midas_db_dir, filetype)
        return _fetch_file_from_s3((s3_file, dest_file))


def _fetch_file_from_s3(packed_args):
    s3_path, dest_file = packed_args
    local_dir = os.path.dirname(dest_file)

    if not os.path.isdir(local_dir):
        command(f"mkdir -p {local_dir}")

    if os.path.exists(dest_file):
        return dest_file

    return download_reference(s3_path, local_dir)
