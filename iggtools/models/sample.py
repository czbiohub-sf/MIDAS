import os
from collections import defaultdict
from iggtools.params.schemas import fetch_schema_by_dbtype, samples_pool_schema, species_profile_schema
from iggtools.common.utils import InputStream, OutputStream, select_from_tsv, command, tsprint


# Executable Documentation

# Low level functions: the Target Files
def get_single_layout(sample_name, dbtype=""):
    def per_species(species_id="", contig_id=""):
        return {
            "species_profile":        f"{sample_name}/species/species_profile.tsv",
            "snps_pileup":            f"{sample_name}/snps/output/{species_id}.snps.tsv.lz4",
            "snps_summary":           f"{sample_name}/snps/output/summary.tsv",
            "genes_coverage":         f"{sample_name}/genes/output/{species_id}.genes.tsv.lz4",
            "genes_summary":          f"{sample_name}/genes/output/summary.tsv",

            "outdir":                 f"{sample_name}/{dbtype}/output",
            "tempdir":                f"{sample_name}/{dbtype}/temp",
            "dbsdir":                 f"{sample_name}/dbs",
            "dbs_tempdir":            f"{sample_name}/dbs/temp",

            "species_alignments_m8":  f"{sample_name}/{dbtype}/temp/alignments.m8",
            "snps_repgenomes_bam":    f"{sample_name}/{dbtype}/temp/repgenomes.bam",
            "genes_pangenomes_bam":   f"{sample_name}/{dbtype}/temp/pangenomes.bam",

            "species_dbs_subdir":     f"{sample_name}/dbs/temp/{species_id}",
            "species_output_subdir":  f"{sample_name}/{dbtype}/output/{species_id}",
            "species_temp_subdir":    f"{sample_name}/{dbtype}/temp/{species_id}",

            "contigs_pileup":         f"{sample_name}/snps/temp/{species_id}/snps_{contig_id}.tsv.lz4",
        }
    return per_species

# Q: how about the base_output_dir? should I put it inside the get_single_layout() function?
# A: No. You should design the low level command for higher level application. Then you
# have a separate command to do the higher level.

# structure things well
# you are the architect, you design it -- this feels really good
# excuatable documentation <- this is super good

# want to add a higher level functions: e.g. pool or orchestral


if False:
    assert os.path.exists(midas_outdir), f"Provided MIDAS output {midas_outdir} for {sample_name} in sample_list is invalid"
    assert os.path.exists(self.data_dir), f"Missing MIDAS {dbtype} directiory for {self.data_dir} for {sample_name}"


class Sample: # pylint: disable=too-few-public-methods
    def __init__(self, sample_name, midas_outdir, dbtype=None):
        self.sample_name = sample_name
        self.midas_outdir = midas_outdir

        self.layout = get_single_layout(sample_name, dbtype=dbtype)
        self.outdir = self.get_target_layout("outdir")
        self.tempdir = self.get_target_layout("tempdir")
        self.dbsdir = self.get_target_layout("dbsdir")

    def get_target_layout(self, filename, species_id="", contig_id=""):
        return os.path.join(self.midas_outdir, self.layout(species_id, contig_id)[filename])

    def create_output_dir(self, debug=False):
        tsprint(f"Create output directory for sample {self.sample_name}.")
        command(f"rm -rf {self.outdir}")
        command(f"mkdir -p {self.outdir}")

        if debug and os.path.exists(self.tempdir):
            tsprint(f"Reusing existing temp data in {self.tempdir} according to --debug flag.")
        else:
            tsprint(f"Create temp directory for sample {self.sample_name}.")
            command(f"rm -rf {self.tempdir}")
            command(f"mkdir -p {self.tempdir}")
            tsprint(f"Create database directory for sample {self.sample_name}.")
            command(f"rm -rf {self.dbsdir}")
            command(f"mkdir -p {self.dbsdir}")

    def create_species_subdir(self, species_ids, debug=False, dirtype=None):
        assert dirtype is not None, f"Need to specify which step the species subdir are built"

        for species_id in species_ids:
            species_subdir = self.get_target_layout(f"species_{dirtype}_subdir", species_id)
            if debug and os.path.exists(species_subdir):
                continue
            command(f"rm -rf {species_subdir}")
            command(f"mkdir -p {species_subdir}")

        ## maybe we should created the outdir and tempdir when know what species will we process
        ## how about provide the index, then we don't those species_id subdirectories .... good question


    def load_species_profile(self, dbtype): # I think this is no longer needed
        species_profile_path = self.layout()["species_profile"]
        assert os.path.exists(species_profile_path), f"Sample::load_species_profile:: missing species profile {species_profile_path} for {self.sample_name}"

    def select_species(self, genome_coverage, species_list=""):
        "Return map of species_id to coverage for the species present in the sample."
        schema = fetch_schema_by_dbtype("species")
        profile = defaultdict()
        with InputStream(self.get_target_layout("species_profile")) as stream:
            for record in select_from_tsv(stream, selected_columns=schema, result_structure=dict):
                if species_list and record["species_id"] not in args.species_list.split(","):
                    continue
                if record["coverage"] >= genome_coverage:
                    profile[record["species_id"]] = record["coverage"]
        return profile

    def remove_output_dir(self):
        command(f"rm -rf {self.tempdir}", check=False)
        command(f"rm -rf {self.outdir}", check=False)
