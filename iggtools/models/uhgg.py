#!/usr/bin/env python3
# A model for the UHGG collection of genomes (aka UHGG database).
import os
from collections import defaultdict
from iggtools.params.outputs import genomes as TABLE_OF_CONTENTS
from iggtools.common.utils import select_from_tsv, sorted_dict, InputStream
from iggtools.params import inputs


# The "output" or built DB layout in S3.
# See https://github.com/czbiohub/iggtools/wiki#target-layout-in-s3

def get_uhgg_layout(species_id, component="", genome_id=""):
    return {
        "genomes_toc":                f"genomes.tsv",

        # External storage of collections of genomes
        "raw_genome_file":            f"{inputs.uhgg_genomes}/{species_id}/{genome_id}.{component}",

        "imported_genome_file":       f"cleaned_imports/{species_id}/{genome_id}/{genome_id}.{component}",
        "imported_genome_log":        f"cleaned_imports/{species_id}/{genome_id}/import_uhgg.log",

        "annotation_file":            f"gene_annotations/{species_id}/{genome_id}/{genome_id}.{component}",
        "annotation_log":             f"gene_annotations/{species_id}/{genome_id}/annotate_genome.log",
        "gene_features_log":          f"gene_annotations/{species_id}/{genome_id}/build_gene_features.log",

        "pangenome_file":             f"pangenomes/{species_id}/{component}",
        "pangenome_log":              f"pangenomes/{species_id}/pangenome_build.log",
        "cluster_info_log":           f"pangenomes/{species_id}/cluster_info.log",

        "marker_genes":               f"markers/{inputs.marker_set}/temp/{species_id}/{genome_id}/{genome_id}.{component}",
        "marker_genes_log":           f"markers/{inputs.marker_set}/temp/{species_id}/{genome_id}/build_marker_genes.log",

        "marker_genes_hmm":           f"markers_models/{inputs.marker_set}/marker_genes.hmm",
        "marker_db_hmm_cutoffs":      f"markers_models/{inputs.marker_set}/marker_genes.mapping_cutoffs{component}",

        "marker_db":                  f"markers/{inputs.marker_set}/{inputs.marker_set}.{component}",
        "build_marker_log":           f"markers/{inputs.marker_set}/build_markerdb.log",

        "marker_centroids":           f"markers/{inputs.marker_set}/marker_centroids/{species_id}.{component}",
        "marker_centroids_log":       f"markers/{inputs.marker_set}/marker_centroids.log",

        "cache_gene_chunks":          f"cache/genes/chunks_chunksize.{genome_id}/{species_id}/{component}.json",
    }


def unified_genome_id(genome_id):
    return "UHGG" + genome_id.replace("GUT_GENOME", "")


def destpath(local_path):
    return os.path.join(inputs.igg, f"{local_path}.lz4")


class UHGG:  # pylint: disable=too-few-public-methods

    def __init__(self, table_of_contents_tsv=""):
        if table_of_contents_tsv:
            self.toc_tsv = table_of_contents_tsv # Read in local genomes.tsv
        else:
            self.toc_tsv = TABLE_OF_CONTENTS # Fetch from S3

        self.species, self.representatives, self.genomes = _UHGG_load(self.toc_tsv)

    def fetch_representative_genome_id(self, species_id):
        return self.representatives[species_id]


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
