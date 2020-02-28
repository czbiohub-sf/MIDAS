import os
from iggtools.params.schemas import fetch_schema_by_dbtype, samples_pool_schema, species_profile_schema
from iggtools.common.utils import InputStream, OutputStream, select_from_tsv, command


# Executable Documentation

# Low level functions: the Target Files
def get_single_layout(sample_name):
    def per_species(species_id=""):
        return {
            "species_profile":  f"{sample_name}/species/species_profile.tsv",
            "snps_pileup":      f"{sample_name}/snps/output/{species_id}.snps.tsv.lz4",
            "snps_summary":     f"{sample_name}/snps/output/summary.tsv",
            "genes_coverage":   f"{sample_name}/genes/output/{species_id}.genes.tsv.lz4",
            "genes_summary":    f"{sample_name}/genes/output/summary.tsv",

            "species_alignments_m8":  f"{sample_name}/species/temp/alignments.m8",
            "snps_repgenomes_bam": f"{sample_name}/snps/temp/repgenomes.bam",
            "genes_pangenomes_bam": f"{sample_name}/genes/temp/pangenomes.bam",
            "dbs_dir": f"{sample_name}/dbs"
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
        self.dbtype = dbtype
        self.outdir = f"{midas_outdir}/{sample_name}/{dbtype}/output"
        self.tempdir = f"{self.outdir}/temp"
        self.layout = get_single_layout(sample_name)
        self.dbsdir = f"{self.outdir}/dbs"

    def create_output_dir(self, debug=False):
        command(f"rm -rf {self.outdir}")
        command(f"mkdir -p {self.outdir}")

        if debug and os.path.exists(self.tempdir):
            tsprint(f"Reusing existing temp data in {self.tempdir} according to --debug flag.")
        else:
            command(f"rm -rf {self.tempdir}")
            command(f"mkdir -p {self.tempdir}")
            command(f"rm -rf {self.dbsdir}")
            command(f"mkdir -p {self.dbsdir}")

        ## I think from the midas_run_genes
        # Add bt2_db_dir to the Sample
        # When we mkdir tempdir, it is okay to add the species subdirectories ...
        ## maybe we should created the outdir and tempdir when know what species will we process
        ## how about provide the index, then we don't those species_id subdirectories .... good question
        #sample.bt2_db_dir = f"{sample.tempdir}/dbs"
        # We should add the dbs/species_id

    def load_species_profile(self, dbtype):
        species_profile_path = self.layout()["species_profile"]
        assert os.path.exists(species_profile_path), f"Sample::load_species_profile:: missing species profile {species_profile_path} for {self.sample_name}"

    def select_species(dbtype, genome_coverage, species_list=""):
        "Return map of species_id to coverage for the species present in the sample."
        schema = fetch_schema_by_dbtype(dbtype)
        profile = defaultdict()
        with InputStream(self.layout()["species_profile"]) as stream:
            for aln in select_from_tsv(stream, selected_columns=schema, result_structure=dict):
                if species_list and aln["species_id"] not in args.species_list.split(","):
                    continue
                if aln["species_coverage"] >= genome_coverage:
                    profile[aln["species_id"]] = aln["species_coverage"]
        return profile

    def remove_output_dir(self):
        command(f"rm -rf {self.tempdir}", check=False)
        command(f"rm -rf {self.outputdir}", check=False)
