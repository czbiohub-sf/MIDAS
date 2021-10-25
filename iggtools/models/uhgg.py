#!/usr/bin/env python3
# A model for the UHGG collection of genomes (aka UHGG database).
import os
from collections import defaultdict
from iggtools.params.outputs import genomes as TABLE_OF_CONTENTS
from iggtools.common.utils import select_from_tsv, sorted_dict, InputStream
from iggtools.params import inputs


def get_uhgg_layout(species_id="", component="", genome_id=""):
    return {
        "raw_genome_file":            f"{inputs.uhgg_genomes}/{species_id}/{genome_id}.{component}",
        "imported_genome":            f"cleaned_imports/{species_id}/{genome_id}/{genome_id}.{component}",
        "imported_genome_log":        f"cleaned_imports/{species_id}/{genome_id}/import_uhgg.log",
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

    def fetch_repgenome_id(self, species_id):
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
