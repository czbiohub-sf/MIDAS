#!/usr/bin/env python3
import os
from midas2.params.schemas import fetch_schema_by_dbtype
from midas2.common.utils import InputStream, select_from_tsv, command, tsprint
from midas2.models.species import filter_species


# Executable Documentation
# Low level functions: the Target Files
def get_single_layout(sample_name, dbtype=""):
    def per_species(species_id="", chunk_id=""):
        return {
            "sample_dir":              f"{sample_name}",
            "outdir":                  f"{sample_name}/{dbtype}",
            "output_subdir":           f"{sample_name}/{dbtype}/{species_id}",

            "tempdir":                 f"{sample_name}/temp/{dbtype}",
            "temp_subdir":             f"{sample_name}/temp/{dbtype}/{species_id}",

            "midas_db_dir":            f"midasdb",
            "bt2_indexes_dir":         f"{sample_name}/bt2_indexes/{dbtype}",

            # species workflow output
            "species_summary":         f"{sample_name}/species/species_profile.tsv",
            "markers_summary":         f"{sample_name}/species/markers_profile.tsv",
            "species_alignments_m8":   f"{sample_name}/temp/species/alignments.m8",
            "species_marker_genes":    f"{sample_name}/temp/species/genes_that_are_marker",
            "species_reads":           f"{sample_name}/temp/species/{species_id}/{chunk_id}.ids",

            # snps workflow output
            "snps_summary":            f"{sample_name}/snps/snps_summary.tsv",
            "snps_pileup":             f"{sample_name}/snps/{species_id}.snps.tsv.lz4",
            "snps_chunk_summary":      f"{sample_name}/snps/chunks_summary.tsv",
            "snps_repgenomes_bam":     f"{sample_name}/temp/snps/repgenomes.bam",
            "species_bam":             f"{sample_name}/temp/snps/{species_id}/{species_id}.bam",
            "species_sorted_bam":      f"{sample_name}/temp/snps/{species_id}/{species_id}.sorted.bam",
            "chunk_pileup":            f"{sample_name}/temp/snps/{species_id}/snps_{chunk_id}.tsv.lz4",

            # genes workflow output
            "genes_summary":           f"{sample_name}/genes/genes_summary.tsv",
            "genes_coverage":          f"{sample_name}/genes/{species_id}.genes.tsv.lz4",
            "genes_chunk_summary":     f"{sample_name}/genes/chunks_summary.tsv",
            "genes_pangenomes_bam":    f"{sample_name}/temp/genes/pangenomes.bam",
            "chunk_coverage":          f"{sample_name}/temp/genes/{species_id}/genes_{chunk_id}.tsv.lz4",
            "chunk_genes_are_markers": f"{sample_name}/temp/genes/{species_id}/markers_{chunk_id}.tsv.lz4",
        }
    return per_species


class Sample: # pylint: disable=too-few-public-methods
    def __init__(self, sample_name, midas_outdir, dbtype=None):
        self.sample_name = sample_name
        self.midas_outdir = midas_outdir
        self.layout = get_single_layout(sample_name, dbtype)
        self.profile = None


    def get_target_layout(self, filename, species_id="", chunk_id=""):
        if isinstance(self.layout(species_id, chunk_id)[filename], list):
            local_file_lists = self.layout(species_id, chunk_id)[filename]
            return [os.path.join(self.midas_outdir, fn) for fn in local_file_lists]
        return os.path.join(self.midas_outdir, self.layout(species_id, chunk_id)[filename])


    def create_dirs(self, list_of_dirnames, debug=False, quiet=False):
        for dirname in list_of_dirnames:
            if dirname == "outdir":
                tsprint(f"Create OUTPUT directory for {self.sample_name}.")
            if dirname == "tempdir":
                tsprint(f"Create TEMP directory for {self.sample_name}.")
            if dirname == "dbsdir":
                tsprint(f"Create DBS directory for {self.sample_name}.")
            create_local_dir(self.get_target_layout(dirname), debug, quiet)


    def create_species_subdirs(self, species_ids, dirname, debug=False, quiet=False):
        for species_id in species_ids:
            species_subdir = self.get_target_layout(f"{dirname}_subdir", species_id)
            create_local_dir(species_subdir, debug, quiet)


    def select_species(self, args, species_list=None):
        """ Parse species profile summary and return list of species for SNPs/Genes analysis """
        if species_list is None:
            species_list = []

        species_profile_fp = self.get_target_layout("species_summary")
        assert os.path.exists(species_profile_fp), f"Need run SPECIES flow before SNPS or GENES for {self.sample_name}"

        schema = fetch_schema_by_dbtype("species")
        assert not any([False for _ in args.select_by.split(',') if _ in schema]), f"Provided {args.select_by} is not in the species profile output for {self.sample_name}"

        species_ids = filter_species(species_profile_fp, args.select_by, args.select_threshold, species_list)
        return species_ids


    def load_profile_by_dbtype(self, dbtype):
        """ Load genes/snps summary in memory and used in Pool model """
        summary_path = self.get_target_layout(f"{dbtype}_summary")
        assert os.path.exists(summary_path), f"load_profile_by_dbtype:: missing {summary_path} for {self.sample_name}"

        schema = fetch_schema_by_dbtype(dbtype)
        profile = dict()
        with InputStream(summary_path) as stream:
            for info in select_from_tsv(stream, selected_columns=schema, result_structure=dict):
                profile[info["species_id"]] = info
        self.profile = profile


    def remove_dirs(self, list_of_dirnames):
        for dirname in list_of_dirnames:
            dirpath = self.get_target_layout(dirname)
            command(f"rm -rf {dirpath}", check=False)


def create_local_dir(dirname, debug, quiet=False):
    if debug and os.path.exists(dirname):
        tsprint(f"Use existing {dirname} according to --debug flag.")
    else:
        command(f"rm -rf {dirname}", quiet)
        command(f"mkdir -p {dirname}", quiet)
