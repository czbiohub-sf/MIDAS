# A model for the UHGG collection of genomes (aka UHGG database).
from collections import defaultdict
from iggtools.params.outputs import genomes as TABLE_OF_CONTENTS
from iggtools.common.utils import select_from_tsv, sorted_dict, InputStream, download_reference, multithreading_map, multiprocessing_map
from iggtools.params import inputs, outputs


def _uhgg_layout(species_id="", component="", genome_id=""):
    return {
        # f"{inputs.uhgg_genomes}/{representative_id}/{genome_id}.fna.lz4"
        "raw_genome_file":            f"{inputs.uhgg_genomes}/{species_id}/{genome_id}.fna.lz4",
        "imported_genome_file":       f"{outputs.cleaned_imports}/{species_id}/{genome_id}/{genome_id}.{component}",
        # s3://microbiome-igg/2.0/pangenomes/GUT_GENOMEDDDDDD/{genes.ffn, centroids.ffn, gene_info.txt}.lz4
        "pangenome_file":             f"{outputs.pangenomes}/{species_id}/{component}",
        "marker_genes_file":          f"{outputs.marker_genes}/phyeco.fa{component}.lz4",
        "marker_genes_mapfile":       f"{outputs.marker_genes}/phyeco.map.lz4",
    }
### old codes, improve the readibility
def raw_genome_file(genome_id, representative_id):
    return f"{inputs.uhgg_genomes}/{representative_id}/{genome_id}.fna.lz4"

### old codes, improve the readibility
def imported_genome_file(genome_id, species_id, component):
    return f"{outputs.cleaned_imports}/{species_id}/{genome_id}/{genome_id}.{component}"
def raw_genome_file(genome_id, representative_id):
    return f"{inputs.uhgg_genomes}/{representative_id}/{genome_id}.fna.lz4"


def get_uhgg_layout(filename, species_id="", component="", genome_id=""):
    return _uhgg_layout(species_id, component, genome_id)[filename]


def fetch_marker_genes(dbs_dir):
    s3_marker_files = []
    for ext in ["", ".bwt", ".header", ".sa", ".sequence"]:
        s3_marker_files.append(get_uhgg_layout("marker_genes_file", species_id="", component=ext))
    s3_marker_files.append(get_uhgg_layout("marker_genes_mapfile"))
    return multithreading_map(fetch_file_from_s3, ((s3file, dbs_dir) for s3file in s3_marker_files))


class UHGG:  # pylint: disable=too-few-public-methods

    def __init__(self, table_of_contents_tsv=TABLE_OF_CONTENTS):
        # Default is to fetch the UHGG table of contents from S3.
        # However, if already downloaded, the local TOC may be passed as a param.
        self.toc_tsv = table_of_contents_tsv
        self.species, self.representatives, self.genomes = _UHGG_load(table_of_contents_tsv)

    def fetch_files(self, list_of_species_ids, dbs_tempdir, filetype=None):
        arguments_list = []
        for species_id in list_of_species_ids:
            local_dir = f"{dbs_tempdir}/{species_id}"
            if filetype == "contigs":
                s3_file = get_uhgg_layout("imported_genome_file", species_id, "fna.lz4", self.representatives[species_id])
            if filetype == "centroids":
                s3_file = get_uhgg_layout("pangenome_file", species_id, "centroids.ffn.lz4")
            if filetype == "genes_info":
                s3_file = get_uhgg_layout("pangenome_file", species_id, "gene_info.txt.lz4")
            arguments_list.append((s3_file, local_dir))
        return multiprocessing_map(fetch_file_from_s3, arguments_list, 20)


    def fetch_representative_genome_id(self, species_id):
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


def fetch_file_from_s3(packed_args):
    s3_path, local_path = packed_args
    return download_reference(s3_path, local_path)
