#!/usr/bin/env python3
import os
from iggtools.params.schemas import fetch_schema_by_dbtype, samples_pool_schema, species_profile_schema, format_data
from iggtools.common.utils import InputStream, OutputStream, select_from_tsv, command, tsprint
from iggtools.models.sample import Sample


def get_pool_layout(dbtype=""):
    def per_species(species_id="", chunk_id=""):
        return {
            "outdir":                f"{dbtype}",
            "outdir_by_species":     f"{dbtype}/{species_id}",

            "tempdir":               f"temp/{dbtype}",
            "tempdir_by_species":    f"temp/{dbtype}/{species_id}",

            "midas_iggdb_dir":       f"midas_iggdb",
            "bt2_indexes_dir":       f"bt2_indexes",

            # species
            "species_prevalence":    f"species/species_prevalence.tsv",
            "species_read_counts":   f"species/species_read_counts.tsv",
            "species_coverage":      f"species/species_coverage.tsv",
            "species_rel_abundance": f"species/species_rel_abundance.tsv",

            # snps
            "snps_summary":          f"snps/snps_summary.tsv",
            "snps_info":             f"snps/{species_id}/{species_id}.snps_info.tsv",
            "snps_freq":             f"snps/{species_id}/{species_id}.snps_freqs.tsv",
            "snps_depth":            f"snps/{species_id}/{species_id}.snps_depth.tsv",
            "snps_info_by_chunk":    f"temp/{dbtype}/{species_id}/cid.{chunk_id}_snps_info.tsv",
            "snps_freq_by_chunk":    f"temp/{dbtype}/{species_id}/cid.{chunk_id}_snps_freqs.tsv",
            "snps_depth_by_chunk":   f"temp/{dbtype}/{species_id}/cid.{chunk_id}_snps_depth.tsv",

            # genes
            "genes_summary":         f"genes/genes_summary.tsv",
            "genes_presabs":         f"genes/{species_id}/{species_id}.genes_presabs.tsv",
            "genes_copynum":         f"genes/{species_id}/{species_id}.genes_copynum.tsv",
            "genes_depth":           f"genes/{species_id}/{species_id}.genes_depth.tsv",
        }
    return per_species


class SamplePool: # pylint: disable=too-few-public-methods

    def __init__(self, samples_list, midas_outdir, dbtype=None):
        self.list_of_samples = samples_list
        self.midas_outdir = midas_outdir
        self.layout = get_pool_layout(dbtype)
        self.samples = self.init_samples(dbtype)


    def get_target_layout(self, filename, species_id="", chunk_id=""):
        return os.path.join(self.midas_outdir, self.layout(species_id, chunk_id)[filename])


    def create_dirs(self, list_of_dirnames, debug=False, quiet=False):
        for dirname in list_of_dirnames:
            if dirname == "outdir":
                tsprint(f"Create OUTPUT directory.")
            if dirname == "tempdir":
                tsprint(f"Create TEMP directory.")
            if dirname == "dbsdir":
                tsprint(f"Create DBS directory.")
            _create_dir(self.get_target_layout(dirname), debug, quiet)


    def create_species_subdirs(self, species_ids, dirname, debug=False, quiet=False):
        dir_to_create = self.get_target_layout(dirname)
        for species_id in species_ids:
            _create_dir(f"{dir_to_create}/{species_id}", debug, quiet)


    def select_species(self, dbtype, args):
        schema = fetch_schema_by_dbtype(dbtype)
        species = {}
        for sample in self.samples:
            with InputStream(sample.get_target_layout(f"{dbtype}_summary")) as stream:
                for record in select_from_tsv(stream, selected_columns=schema, result_structure=dict):
                    species_id = record["species_id"]
                    # Skip unspeficied species
                    if (args.species_list and species_id not in args.species_list.split(",")):
                        continue
                    # Read in all the species_id in the species profile
                    if species_id not in species:
                        species[species_id] = Species(species_id)
                    # Skip low-coverage <species, sample>
                    if record['mean_coverage'] < args.genome_depth:
                        continue
                    # Skip low prevalent <species, sample>
                    if (dbtype == "snps" and record['fraction_covered'] < args.genome_coverage):
                        continue
                    # Select high quality sample-species pairs
                    species[species_id].samples.append(sample)

        list_of_species = list(species.values())
        # Sort list_of_species by samples_count in descending order
        list_of_species = _sort_species(list_of_species)
        # Second round of filters based on prevalence (sample_counts)
        list_of_species = _filter_species(list_of_species, args)
        return {species.id:species for species in list_of_species}


    def init_samples(self, dbtype):
        """ read in table-of-content: samples_list """
        samples = []
        with InputStream(self.list_of_samples) as stream:
            for row in select_from_tsv(stream, selected_columns=samples_pool_schema, result_structure=dict):
                sample = Sample(row["sample_name"], row["midas_outdir"], dbtype)
                # load profile_summary into memory for easy access
                sample.load_profile_by_dbtype(dbtype)
                samples.append(sample)
        return samples


    def fetch_samples_names(self):
        return [sample.sample_name for sample in self.samples]


    def write_summary_files(self, dict_of_species, dbtype):
        """ Write snps/genes summary files for current samples pool """

        summary_file = self.get_target_layout(f"{dbtype}_summary")
        summary_header = list(fetch_schema_by_dbtype(dbtype).keys())[1:]
        with OutputStream(summary_file) as stream:
            stream.write("\t".join(["species_id", "sample_name"] + summary_header) + "\n")
            for species in dict_of_species.values():
                for sample in species.samples:
                    row = list(sample.profile[species.id].values())
                    stream.write("\t".join([str(row[0]), sample.sample_name]) + "\t" + "\t".join(map(format_data, row[1:])) + "\n")


    def remove_dirs(self, list_of_dirnames):
        for dirname in list_of_dirnames:
            dirpath = self.get_target_layout(dirname)
            command(f"rm -rf {dirpath}", check=False)


class Species:
    """ Base class for species """
    def __init__(self, id):
        self.id = id
        self.samples = []
        self.samples_depth = []


    def fetch_samples_depth(self):
        return [sample.profile[self.id]["mean_coverage"] for sample in self.samples]


    def fetch_samples_names(self):
        list_of_sample_objects = list(self.samples)
        return [sample.sample_name for sample in list_of_sample_objects]


def _sort_species(species):
    """ Sort list_of_species by samples_count in descending order """
    species_sorted = sorted(((sp, len(sp.samples)) for sp in species), key=lambda x: x[1], reverse=True)
    return [sp[0] for sp in species_sorted]


def _filter_species(species, args):
    """ Filter out low prevalent species using samples_count cutoff """
    species_keep = []
    for sp in species:
        sp.samples_count = len(sp.samples)
        # skip low prevalent species
        if sp.samples_count < args.sample_counts:
            continue
        sp.samples_depth = sp.fetch_samples_depth()
        species_keep.append(sp)
    return species_keep


def _create_dir(dirname, debug, quiet=False):
    if debug and os.path.exists(dirname):
        tsprint(f"Use existing {dirname} according to --debug flag.")
    else:
        command(f"rm -rf {dirname}", quiet)
        command(f"mkdir -p {dirname}", quiet)
