import os
from collections import defaultdict
from iggtools.params.schemas import fetch_schema_by_dbtype, samples_pool_schema, species_profile_schema
from iggtools.common.utils import InputStream, OutputStream, select_from_tsv, command
from iggtools.models.sample import Sample


def get_pool_layout(dbtype=""):
    def per_species(species_id=""):
        return {
            "species_prevalence":    f"species/species_prevalence.tsv",

            "snps_summary":          f"snps/output/snps_summary.tsv",
            "snps_info":             f"snps/output/{species_id}/{species_id}.snps_info.tsv",
            "snps_freq":             f"snps/output/{species_id}/{species_id}.snps_freqs.tsv",
            "snps_depth":            f"snps/output/{species_id}/{species_id}.snps_depth.tsv",

            "genes_summary":         f"genes/output/summary.tsv",
            "genes_presabs":         f"genes/output/{species_id}/{species_id}.genes_presabs.tsv",
            "genes_copynum":         f"genes/output/{species_id}/{species_id}.genes_copynum.tsv",
            "genes_depth":           f"genes/output/{species_id}/{species_id}.genes_depth.tsv",

            "outdir":                f"{dbtype}/output",
            "tempdir":               f"{dbtype}/temp",
            "dbsdir":                f"{dbtype}/temp/dbs",
        }
    return per_species


class Pool: # pylint: disable=too-few-public-methods

    def __init__(self, samples_list, midas_outdir, dbtype=None):
        self.list_of_samples = samples_list
        self.midas_outdir = midas_outdir

        self.samples = self.init_samples(dbtype)
        self.layout = get_pool_layout(dbtype=dbtype)
        self.outdir = self.get_target_layout("outdir")
        self.tempdir = self.get_target_layout("tempdir")
        self.dbsdir = self.get_target_layout("dbsdir")


    def get_target_layout(self, filename, species_id=""):
        return os.path.join(self.midas_outdir, self.layout(species_id)[filename])

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

    def create_species_subdir(self, species_ids, debug=False):
        dbs_tempdir = self.get_target_layout("dbs_tempdir")
        for species_id in species_ids:
            if debug and os.path.exists(f"{dbs_tempdir}/{species_id}"):
                continue
            command(f"rm -rf {dbs_tempdir}/{species_id}")
            command(f"mkdir -p {dbs_tempdir}/{species_id}")


    def init_samples(self, dbtype):
        """ read in table-of-content: sample_name\tpath/to/midas_output """
        samples = []
        with InputStream(self.list_of_samples) as stream:
            for row in select_from_tsv(stream, selected_columns=samples_pool_schema, result_structure=dict):
                sample = Sample(row["sample_name"], row["midas_outdir"], dbtype)
                samples.append(sample)
        return samples

    def create_outdir(self, outdir, dbtype):
        # only created once in the beginning
        outdir = f"{outdir}/merged/{dbtype}"
        command(f"rm -rf {outdir}")
        command(f"mkdir -p {outdir}")
        self.outdir = outdir

    def create_tempdir(self, paramstr):
        # only created once in the beginning
        tempdir = f"{self.outdir}/temp_{paramstr}"
        command(f"rm -rf {tempdir}")
        command(f"mkdir -p {tempdir}")
        self.tempdir = tempdir

    def create_species_tempdir(self, species_id):
        # first created when creating lookup tables
        species_tempdir = f"{self.tempdir}/{species_id}"
        if not os.path.exists(species_tempdir):
            command(f"mkdir -p {species_tempdir}")
        return species_tempdir

    def create_species_lookup_table(self, species_id):
        species_tempdir = self.create_species_tempdir(species_id)
        proc_lookup_file = f"{species_tempdir}/procid_lookup.tsv"
        return proc_lookup_file

    def create_species_pool_snps_chunk(self, species_id, process_id):
        # chunk file for one pair of species_id-process_id
        species_tempdir = self.create_species_tempdir(species_id)
        snps_info_fp = f"{species_tempdir}/proc.{process_id}_snps_info.tsv"
        snps_freq_fp = f"{species_tempdir}/proc.{process_id}_snps_freqs.tsv"
        snps_depth_fp = f"{species_tempdir}/proc.{process_id}_snps_depth.tsv"
        return (snps_info_fp, snps_freq_fp, snps_depth_fp)

    def create_species_outdir(self, species_id):
        # first created when merge_chunks_by_species
        species_outdir = f"{self.outdir}/{species_id}"
        if not os.path.exists(species_outdir):
            command(f"mkdir -p {species_outdir}")
        return species_outdir

    def create_species_pool_snps_results(self, species_id):
        species_outdir = self.create_species_outdir(species_id)
        snps_info_fp = f"{species_outdir}/{species_id}_snps_info.tsv"
        snps_freq_fp = f"{species_outdir}/{species_id}_freqs.tsv"
        snps_depth_fp = f"{species_outdir}/{species_id}_snps_depth.tsv"
        return (snps_info_fp, snps_freq_fp, snps_depth_fp)

    def fetch_sample_names(self):
        sample_names = [sample.sample_name for sample in self.samples]
        return sample_names

    def fetch_merged_species_output_files(self):
        cols_to_write = list(species_profile_schema.keys())[1:]
        return {col: f"{self.outdir}/{col}.tsv" for col in cols_to_write}

    def fetch_species_prevalence_path(self):
        return f"{self.outdir}/species_prevalence.tsv"


class Species:
    """ Base class for species """
    def __init__(self, id):
        self.id = id
        self.samples = []
        self.samples_depth = []

    def fetch_samples_depth(self):
        return [sample.profile[self.id]["mean_coverage"] for sample in self.samples]

    def write_summary(self, dbtype, outdir):
        """ Write snps/genes summary files for current samples pool """
        summary_dir = f"{outdir}/{self.id}/{dbtype}_summary.tsv"
        with OutputStream(summary_dir) as stream:
            header = ["species_id", "sample_name"] + list(fetch_schema_by_dbtype(dbtype).keys())[1:]
            stream.write("\t".join(header) + "\n")
            for sample in self.samples:
                record = [self.id, sample.sample_name] + list(sample.profile[self.id].values())[1:]
                stream.write("\t".join(map(str, record)) + "\n")

    def fetch_samples_name(self):
        list_of_sample_objects = list(self.samples)
        samples_names = [sample.sample_name for sample in list_of_sample_objects]
        return samples_names


def init_species(pool_of_samples, dbtype, args):
    schema = fetch_schema_by_dbtype(dbtype)
    species = {}
    for sample in pool_of_samples.samples:
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

    return list(species.values())


def sort_species(species):
    """ Sort list_of_Species by samples_count in descending order """
    species_sorted = sorted(((sp, len(sp.samples)) for sp in species), key=lambda x: x[1], reverse=True)
    return [sp[0] for sp in species_sorted]


def filter_species(species, args):
    """ Filter out low prevalent species using samples_count cutoff """
    species_keep = []
    for sp in species:
        sp.samples_count = len(sp.samples)
        if sp.samples_count < args.sample_counts:
            continue # skip low prevalent species
        sp.samples_depth = sp.fetch_samples_depth()
        species_keep.append(sp)
    return species_keep


def select_species(pool_of_samples, dbtype, args):
    species = init_species(pool_of_samples, dbtype, args)
    species = sort_species(species)
    species = filter_species(species, args)
    return species


def search_species(list_of_species, species_id):
    curr_species = [species for species in list_of_species if species.id == species_id]
    assert len(curr_species) == 1, f"pool::search_species duplicated Species in list_of_species"
    return curr_species[0]
