from iggtools.models.sample import Sample
from iggtools.params.schemas import fetch_schema_by_dbtype, samples_pool_schema
from iggtools.common.utils import InputStream, OutputStream, select_from_tsv


class Pool: # pylint: disable=too-few-public-methods

    def __init__(self, sample_list, dbtype=None):
        self.pool_tsv = sample_list
        self.samples = self.init_samples(dbtype)

    def init_samples(self, dbtype):
        """ read in table-of-content: sample_name\tpath/to/midas_output """
        samples = []
        with InputStream(self.pool_tsv) as stream:
            for row in select_from_tsv(stream, selected_columns=samples_pool_schema, result_structure=dict):
                sample = Sample(row["sample_name"], row["midas_outdir"], dbtype)
                samples.append(sample)
        return samples


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


def _filter_sample_species(info, args, dbtype):
    """ select high quality sample-species pairs"""
    if info['mean_coverage'] < args.genome_depth:
        return False # skip low-coverage <species, sample>
    if (args.species_id and info["species_id"] not in args.species_id.split(",")):
        return False # skip unspeficied species
    if (dbtype == "snps" and info['fraction_covered'] < args.genome_coverage):
        return False # skip low prevalent <species, sample>
    return True


def init_species(pool_of_samples, dbtype, args):
    species = {}
    for sample in pool_of_samples.samples:
        for info in sample.profile.values():
            species_id = info["species_id"]
            if species_id not in species:
                species[species_id] = Species(species_id)
            if _filter_sample_species(info, args, dbtype):
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


def select_species(pool_of_samples, args, dbtype):
    species = init_species(pool_of_samples, dbtype, args)
    species = sort_species(species)
    species = filter_species(species, args)
    return species
