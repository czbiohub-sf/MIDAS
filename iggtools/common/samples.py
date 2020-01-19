#!/usr/bin/env python3
from iggtools.common.utils import InputStream, select_from_tsv


def parse_species_profile(outdir):
    "Return map of species_id to coverage for the species present in the sample."
    with InputStream(f"{outdir}/species/species_profile.txt") as stream:
        return dict(select_from_tsv(stream, {"species_id": str, "coverage": float}))


def select_species(species_profile, coverage_threshold):
    return {species_id: species_coverage for species_id, species_coverage in species_profile.items() if species_coverage >= coverage_threshold}
