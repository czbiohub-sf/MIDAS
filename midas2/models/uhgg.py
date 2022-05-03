#!/usr/bin/env python3
from collections import defaultdict
from midas2.params.outputs import genomes as TABLE_OF_CONTENTS
from midas2.common.utils import select_from_tsv, sorted_dict, InputStream, OutputStream, find_files, tsprint
from midas2.params import inputs, outputs


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
