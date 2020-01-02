# A model for the UHGG collection of genomes (aka UHGG database).
from collections import defaultdict
from iggtools.params.outputs import genomes as TABLE_OF_CONTENTS
from iggtools.common.utils import select_from_tsv, sorted_dict, InputStream


class UHGG:  # pylint: disable=too-few-public-methods

    def __init__(self, table_of_contents_tsv=TABLE_OF_CONTENTS):
        # Default is to fetch the UHGG table of contents from S3.
        # However, if already downloaded, the local TOC may be passed as a param.
        self.toc_tsv = table_of_contents_tsv
        self.species, self.representatives, self.genomes = _UHGG_load(table_of_contents_tsv)


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
