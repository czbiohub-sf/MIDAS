#!/usr/bin/env python3
# A model for the UHGG collection of genomes (aka UHGG database).
from collections import defaultdict
from midas2.params.outputs import genomes as TABLE_OF_CONTENTS
from midas2.common.utils import select_from_tsv, sorted_dict, InputStream, OutputStream, find_files, tsprint
from midas2.params import inputs, outputs


def get_uhgg_layout(species_id="", component="", genome_id=""):
    return {
        "raw_genome_file":            f"{inputs.uhgg_genomes}/{species_id}/{genome_id}.{component}",
        "imported_genome":            f"cleaned_imports/{species_id}/{genome_id}/{genome_id}.{component}",
        "imported_genome_log":        f"cleaned_imports/{species_id}/{genome_id}/import_uhgg.log",
    }


def unified_genome_id(genome_id):
    # this is an obsolete tag only for UHGG.
    return "UHGG" + genome_id.replace("GUT_GENOME", "")


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


def build_uhgg_toc(args):
    # Moved from init.init
    toc_genome = outputs.genomes("uhgg")

    msg = f"Building {toc_genome}."
    if find_files(toc_genome):
        if not args.force:
            tsprint(f"Destination {toc_genome} already exists.  Specify --force to overwrite.")
            return
        msg = f"Rebuilding {toc_genome}."
    tsprint(msg)

    id_remap = {}
    with InputStream(inputs.alt_species_ids) as ids:
        for row in select_from_tsv(ids, selected_columns=["alt_species_id", "species_id"]):
            new_id, old_id = row
            id_remap[old_id] = new_id

    seen_genomes, seen_species = set(), set()
    with OutputStream(toc_genome) as out:

        target_columns = ["genome", "species", "representative", "genome_is_representative"]
        out.write("\t".join(target_columns) + "\n")

        with InputStream(inputs.genomes2species) as g2s:
            for row in select_from_tsv(g2s, selected_columns=["MAG_code", "Species_id"]):
                genome, representative = row
                species = id_remap[representative]
                genome_is_representative = str(int(genome == representative))
                target_row = [genome, species, representative, genome_is_representative]
                out.write("\t".join(target_row) + "\n")
                seen_genomes.add(genome)
                seen_species.add(species)

    tsprint(f"Emitted {len(seen_genomes)} genomes and {len(seen_species)} species to {outputs.genomes}.")
