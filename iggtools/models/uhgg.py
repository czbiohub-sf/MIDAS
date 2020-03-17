# A model for the UHGG collection of genomes (aka UHGG database).
import os
from collections import defaultdict
from iggtools.params.outputs import genomes as TABLE_OF_CONTENTS
from iggtools.common.utils import select_from_tsv, sorted_dict, InputStream, download_reference, multithreading_map
from iggtools.params import inputs, outputs


def get_uhgg_layout(species_id="", component=""):
    return {
        "marker_genes_mapfile": f"{outputs.marker_genes}/phyeco.map.lz4",
        # s3://microbiome-igg/2.0/pangenomes/GUT_GENOMEDDDDDD/{genes.ffn, centroids.ffn, gene_info.txt}
        "pangenome_file": f"{outputs.pangenomes}/{species_id}/{component}",
    }

class UHGG:  # pylint: disable=too-few-public-methods

    def __init__(self, table_of_contents_tsv=TABLE_OF_CONTENTS):
        # Default is to fetch the UHGG table of contents from S3.
        # However, if already downloaded, the local TOC may be passed as a param.
        self.toc_tsv = table_of_contents_tsv
        self.species, self.representatives, self.genomes = _UHGG_load(table_of_contents_tsv)

    def fetch_representative_genome_id(self, species_id):
        return self.representatives[species_id]

    def fetch_marker_genes(self, dbs_dir):
        target_markers_db_files = [
            f"{outputs.marker_genes}/phyeco.fa{ext}.lz4" for ext in ["", ".bwt", ".header", ".sa", ".sequence"]
            ] + \
            [f"{outputs.marker_genes}/phyeco.map.lz4"]
        argument_list = ((mfile, dbs_dir) for mfile in target_markers_db_files)
        markers_db_files = multithreading_map(fetch_file_from_s3, argument_list)
        return markers_db_files

    def fetch_contigs(self, species_ids, bt2_db_dir):
        argument_list = []
        for species_id in species_ids:
            assert os.path.exists(f"{bt2_db_dir}/{species_id}"), "Fail to create {bt2_db_dir}/{species_id} in create_species_subdir()"
            argument_list.append((imported_genome_file(self.representatives[species_id], species_id, "fna.lz4"), f"{bt2_db_dir}/{species_id}"))
        contigs_files = multithreading_map(fetch_file_from_s3, argument_list, num_threads=20)
        return contigs_files

    def fetch_centroids(self, species_ids, dbs_tempdir):
        argument_list = []
        for species_id in species_ids:
            centroid_file = get_uhgg_layout(species_id, "centroids.ffn.lz4")["pangenome_file"]
            print(centroid_file)
            argument_list.append((centroid_file, f"{dbs_tempdir}/{species_id}"))
        centroids_files = multithreading_map(fetch_file_from_s3, argument_list, num_threads=20)
        return centroids_files

    def fetch_genes_info(self, species_ids, dbs_tempdir):
        argument_list = []
        for species_id in species_ids:
            genes_info_file = get_uhgg_layout(species_id, "gene_info.txt.lz4")["pangenome_file"]
            argument_list.append((genes_info_file, f"{dbs_tempdir}/{species_id}"))
        genes_info_files = multithreading_map(fetch_file_from_s3, argument_list, num_threads=20)
        return genes_info_files


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


def imported_genome_file(genome_id, species_id, component):
    return f"{outputs.cleaned_imports}/{species_id}/{genome_id}/{genome_id}.{component}"

def raw_genome_file(genome_id, representative_id):
    return f"{inputs.uhgg_genomes}/{representative_id}/{genome_id}.fna.lz4"

def pangenome_file(species_id, component):
    # s3://microbiome-igg/2.0/pangenomes/GUT_GENOMEDDDDDD/{genes.ffn, centroids.ffn, gene_info.txt}
    return f"{outputs.pangenomes}/{species_id}/{component}"

def fetch_file_from_s3(packed_args):
    s3_path, local_path = packed_args
    return download_reference(s3_path, local_path)
