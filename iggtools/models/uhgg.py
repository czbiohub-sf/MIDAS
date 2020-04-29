# A model for the UHGG collection of genomes (aka UHGG database).
import os
from collections import defaultdict
from iggtools.params.outputs import genomes as TABLE_OF_CONTENTS
from iggtools.common.utils import select_from_tsv, sorted_dict, InputStream, download_reference, multiprocessing_map, multithreading_map, command
from iggtools.params import outputs, inputs
from iggtools.params.inputs import igg

MARKER_FILE_EXTS = ["fa", "fa.bwt", "fa.header", "fa.sa", "fa.sequence", "map"]


def get_uhgg_layout(filename, species_id="", component="", genome_id=""):
    return {
        "genomes_toc":                f"genomes.tsv",
        # f"{inputs.uhgg_genomes}/{representative_id}/{genome_id}.fna.lz4"
        #"raw_genome_file":            f"{inputs.uhgg_genomes}/{species_id}/{genome_id}.fna.lz4",
        #"imported_genome_file":       f"{outputs.cleaned_imports}/{species_id}/{genome_id}/{genome_id}.{component}",
        "imported_genome_file":       f"cleaned_imports/{species_id}/{genome_id}/{genome_id}.{component}",
        # 100001/{genes.ffn, centroids.ffn, gene_info.txt}.lz4
        "pangenome_file":             f"pangenomes/{species_id}/{component}",

        "marker_db":                  f"marker_genes/phyeco/phyeco.{component}.lz4",
    }

### old codes, improve the readibility
## there two functions were used in import_uhgg.py leave it alone for now
def imported_genome_file(genome_id, species_id, component):
    return f"{outputs.cleaned_imports}/{species_id}/{genome_id}/{genome_id}.{component}"
def raw_genome_file(genome_id, representative_id):
    return f"{inputs.uhgg_genomes}/{representative_id}/{genome_id}.fna.lz4"


def fetch_marker_genes(dbs_dir):
    s3_marker_files = []
    for ext in MARKER_FILE_EXTS:
        s3_marker_files.append(get_uhgg_layout("marker_genes_file", species_id="", component=ext))
    s3_marker_files.append(get_uhgg_layout("marker_genes_mapfile"))
    return multithreading_map(fetch_file_from_s3, ((s3file, dbs_dir) for s3file in s3_marker_files))


class MIDAS_IGGDB: # pylint: disable=too-few-public-methods

    def __init__(self, midas_iggdb_dir):
        self.local_toc = download_reference(outputs.genomes, midas_iggdb_dir)
        self.midas_iggdb_dir = midas_iggdb_dir
        self.uhgg = UHGG(self.local_toc)

    def get_target_layout(self, filename, species_id="", component="", genome_id="", remote=False):
        file_name = get_uhgg_layout(filename, species_id, component, genome_id)
        target_filepath = f"{igg}/{file_name}" if remote else os.path.join(self.midas_iggdb_dir, file_name)
        return target_filepath


    def fetch_files(self, list_of_species_ids, filetype=None):
        # Fetch igg/2.0 files to local midas_iggdb


        for species_id in list_of_species_ids:

            args_dict = {}
            fetched_files = {}

            if filetype == "marker_db":
                for ext in MARKER_FILE_EXTS:
                    s3_file = self.get_target_layout("marker_db", species_id="", component=ext, genome_id="", True)
                    local_file = self.get_target_layout("marker_db", species_id="", component=ext, genome_id="", False)
                    fetch_files[ext] = fetch_file_from_s3(s3_file, local_file)
                return fetched_files

            if filetype == "contigs":
                s3_file = self.get_target_layout("imported_genome_file", species_id, "fna.lz4", self.uhgg.representatives[species_id], True)
                local_file = self.get_target_layout("imported_genome_file", species_id, "fna", self.uhgg.representatives[species_id], False)

            if filetype == "centroids":
                s3_file = self.get_target_layout("pangenome_file", species_id, "centroids.ffn.lz4", True)
                local_file = self.get_target_layout("pangenome_file", species_id, "centroids.ffn", True)

            if filetype == "genes_info":
                s3_file = self.get_target_layout("pangenome_file", species_id, "gene_info.txt.lz4", True)
                local_file = self.get_target_layout("pangenome_file", species_id, "gene_info.txt", False)

            if not os.path.exists(local_file):
                local_dir = os.path.dirname(local_file)
                command(f"mkdir -p {local_dir}")
                args_dict[species_id] = ((s3_file, local_dir))
            else:
                fetched_files[species_id] = local_file

        _fetched_files = multiprocessing_map(fetch_file_from_s3, list(args_dict.values()), 20)
        for species_index, species_id in enumerate(args_dict.keys()):
            fetched_files[species_id] = _fetched_files[species_index]
        return fetched_files


class UHGG:  # pylint: disable=too-few-public-methods

    def __init__(self, table_of_contents_tsv=TABLE_OF_CONTENTS):
        # Default is to fetch the UHGG table of contents from S3.
        # However, if already downloaded, the local TOC may be passed as a param.
        self.toc_tsv = table_of_contents_tsv
        self.species, self.representatives, self.genomes = _UHGG_load(table_of_contents_tsv)


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
