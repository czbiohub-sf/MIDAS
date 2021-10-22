#!/usr/bin/env python3
# A model for the UHGG collection of genomes (aka UHGG database).
import os
from iggtools.common.utils import download_reference, command, multiprocessing_map, tsprint, InputStream, select_from_tsv, drop_lz4
from iggtools.models.uhgg import UHGG, get_uhgg_layout, destpath
from iggtools.params.inputs import MARKER_FILE_EXTS, igg_dict, marker_set


# The "output" or built DB layout in S3.
# See https://github.com/czbiohub/iggtools/wiki#target-layout-in-s3
def get_midasdb_layout(species_id="", genome_id="", component=""):
    return {
        "table_of_contents":             f"genomes.tsv",

        "imported_genome":               f"cleaned_imports/{species_id}/{genome_id}/{genome_id}.{component}",
        "imported_genome_log":           f"cleaned_imports/{species_id}/{genome_id}/import_uhgg.log",

        "representative_genome":         f"gene_annotations/{species_id}/{genome_id}/{genome_id}.fna",

        "marker_db":                     [f"markers/{marker_set}/{marker_set}.{ext}" for ext in MARKER_FILE_EXTS],
        "marker_db_hmm":                 f"markers_models/{marker_set}/marker_genes.hmm",
        "marker_db_hmm_cutoffs":         f"markers_models/{marker_set}/marker_genes.mapping_cutoffs",
        "build_markerdb_log":            f"markers/{marker_set}/build_markerdb.log",

        "pangenome_file":                f"pangenomes/{species_id}/{component}",
        "pangenome_centroids":           f"pangenomes/{species_id}/centroids.ffn",
        "pangenome_genes_info":          f"pangenomes/{species_id}/gene_info.txt",
        "pangenome_genes_len":           f"pangenomes/{species_id}/genes.len",
        "pangenome_cluster_info":        f"pangenomes/{species_id}/cluster_info.txt",
        "cluster_info_log":              f"pangenomes/{species_id}/cluster_info.log",

        "chunks_centroids":              f"chunks/genes/chunksize.{component}/{species_id}.json",
        "chunks_sites_run":              f"chunks/sites/run/chunksize.{component}/{species_id}/{genome_id}.json",
        "chunks_sites_merge":            f"chunks/sites/merge/chunksize.{component}/{species_id}/{genome_id}.json",

        "marker_genes":                  f"markers/{marker_set}/temp/{species_id}/{genome_id}/{genome_id}.{component}",
        "marker_genes_seq":              f"markers/{marker_set}/temp/{species_id}/{genome_id}/{genome_id}.markers.fa",
        "marker_genes_map":              f"markers/{marker_set}/temp/{species_id}/{genome_id}/{genome_id}.markers.map",
        "marker_genes_hmmsearch":        f"markers/{marker_set}/temp/{species_id}/{genome_id}/{genome_id}.hmmsearch",
        "marker_genes_log":              f"markers/{marker_set}/temp/{species_id}/{genome_id}/build_marker_genes.log",

        "annotation_file":               f"gene_annotations/{species_id}/{genome_id}/{genome_id}.{component}",
        "annotation_faa":                f"gene_annotations/{species_id}/{genome_id}/{genome_id}.faa",
        "annotation_ffn":                f"gene_annotations/{species_id}/{genome_id}/{genome_id}.ffn",
        "annotation_fna":                f"gene_annotations/{species_id}/{genome_id}/{genome_id}.fna",
        "annotation_gff":                f"gene_annotations/{species_id}/{genome_id}/{genome_id}.gff",
        "annotation_tsv":                f"gene_annotations/{species_id}/{genome_id}/{genome_id}.tsv",
        "annotation_log":                f"gene_annotations/{species_id}/{genome_id}/annotate_genome.log",
        "annotation_genes":              f"gene_annotations/{species_id}/{genome_id}/{genome_id}.genes",

        "gene_features_log":             f"gene_annotations/{species_id}/{genome_id}/build_gene_features.log",
    }


class MIDAS_DB: # pylint: disable=too-few-public-methods

    def __init__(self, db_dir=".", db_name="uhgg", num_cores=1):
        self.db_dir = os.path.abspath(db_dir) # local dir
        self.db_name = igg_dict[db_name] # s3 dir
        self.num_cores = num_cores
        self.local_toc = self.fetch_files("table_of_contents")

        self.uhgg = UHGG(self.local_toc)


    def get_repgenome_id(self, species_id):
        return self.uhgg.representatives[species_id]


    def construct_local_path(self, filename, species_id="", genome_id="", component=""):
        local_path = _get_local_path(get_midasdb_layout(species_id, genome_id, component)[filename], self.db_dir)
        return local_path


    def construct_dest_path(self, filename, species_id="", genome_id="", component=""):
        s3_path = _get_dest_path(get_midasdb_layout(species_id, genome_id, component)[filename], self.db_name)
        return s3_path


    def construct_file_tuple(self, filename, species_id="", genome_id="", component=""):
        local_path = self.construct_local_path(filename, species_id, genome_id, component)
        s3_path = self.construct_dest_path(filename, species_id, genome_id, component)
        return (s3_path, local_path)


    def fetch_file(self, filename, species_id="", genome_id="", component=""):
        (s3_path, local_path) = self.construct_file_tuple(filename, species_id, genome_id, component)
        return _fetch_file_from_s3((s3_path, local_path))


    def fetch_files(self, filename, list_of_species=""):
        """ Fetch files from S3 to local MIDAS_DB """
        if list_of_species:
            args_list = []
            for species_id in list_of_species:
                genome_id = self.get_repgenome_id(species_id)
                args_list.append(self.construct_file_tuple(filename, species_id, genome_id))

            if len(list_of_species) > 1:
                _fetched_files = multiprocessing_map(_fetch_file_from_s3, args_list, self.num_cores) #<-----
            else:
                _fetched_files = [_fetch_file_from_s3(args_list[0])]

            fetched_files = {}
            for species_index, species_id in enumerate(list_of_species):
                fetched_files[species_id] = _fetched_files[species_index]
            return fetched_files

        if filename == "marker_db":
            l1 = self.construct_dest_path(filename)
            l2 = self.construct_local_path(filename)
            list_of_file_tuples = list(zip(l1, l2))
            fetched_files = [_fetch_file_from_s3(_) for _ in list_of_file_tuples]
            return dict(zip(MARKER_FILE_EXTS, fetched_files))

        # Fetch single file with filename as the uhgg layout key
        if filename in get_uhgg_layout(species_id=""):
            s3_path = self.get_target_layout(filename, True)
            local_path = self.get_target_layout(filename, False)
            #s3_path = self.construct_dest_path(filename, True)
            #local_path = self.construct_local_path(filename, False)
        else:
            (s3_path, local_path) = self.construct_file_tuple(filename)
        return _fetch_file_from_s3((s3_path, local_path))


    def get_target_layout_old(self, filename, remote, component="", species_id="", genome_id=""):
        ## TODO: (1) left over outside this function (2) keys in uhgg_layout but not in midasdb_layout
        local_path = get_uhgg_layout(species_id, component, genome_id)[filename]
        target_path = destpath(local_path) if remote else os.path.join(self.db_dir, local_path)
        return target_path


    def get_target_layout(self, filename, remote, species_id="", genome_id="", component=""):
        local_path = self.construct_local_path(filename, species_id, genome_id, component)
        s3_path = self.construct_dest_path(filename, species_id, genome_id, component)
        target_path = s3_path if remote else local_path
        return target_path


    def test_markers(self):
        # this is only the structured file, no absoluate path
        species_id = "117086" #117082
        genome_id = "GCA_900552055.1" #GCA_003451595.1
        chunk_size = 50000

        list_of_genomes = self.uhgg.species[species_id].keys()
        list_of_mapfiles = [self.fetch_file("marker_genes_map", species_id, genome_id) for genome_id in list_of_genomes]
        #cat_files(list_of_mapfiles, "mapfile", 20)

        tsprint("1")
        filter_cmd = f"awk \'$8 != \"\"\'"
        cluster_info_fp = self.fetch_file("pangenome_cluster_info", species_id)
        print(cluster_info_fp)
        list_of_marker_genes = []
        dict_of_genes_are_markers = dict()
        with InputStream(cluster_info_fp, filter_cmd) as stream:
            for r in select_from_tsv(stream, selected_columns=["centroid_99", "marker_id"], result_structure=dict):
                dict_of_genes_are_markers[r["centroid_99"]] = r # r only has two elements: centroids_99 and marker_id
                if r["marker_id"] not in list_of_marker_genes:
                    list_of_marker_genes.append(r["marker_id"])
        print(dict_of_genes_are_markers)
        tsprint("2")

    def test(self):
        species_id = "100001" #117082
        genome_id = "GUT_GENOME000001" #GCA_003451595.1
        print(destpath(get_uhgg_layout(species_id, "fna", genome_id)["imported_genome"]))
        print(self.get_target_layout("imported_genome", True, species_id, genome_id, "fna"))
        print(self.fetch_file("imported_genome", species_id, genome_id, "fna"))

        dest_file = self.get_target_layout("annotation_fna", True, species_id, genome_id)
        print(self.get_target_layout("annotation_file", True, species_id, genome_id, "fna"))
        last_output = drop_lz4(os.path.basename(dest_file))
        print(last_output)
        print(dest_file)

        output_files = [
            f"{genome_id}.faa",
            f"{genome_id}.ffn",
            f"{genome_id}.fna",
            f"{genome_id}.gff",
            f"{genome_id}.tsv"
        ]
        upload_tasks = []
        upload_tasks2 = []
        for o in output_files:
            #olz = o + ".lz4"
            otype = o.rsplit(".")[-1]
            if o != last_output:
                upload_tasks2.append((o, self.get_target_layout("annotation_file", True, species_id, genome_id, otype)))
                upload_tasks.append((o, destpath(get_uhgg_layout(species_id, otype, genome_id)["annotation_file"])))
        print(upload_tasks)
        print("\n\n")
        print(upload_tasks2)


def _get_dest_path(file_name, db_name):
    """ Append S3 MIDAS_DB path to file(s) """
    if isinstance(file_name, list):
        return [os.path.join(db_name, f"{fn}.lz4") for fn in file_name]
    return os.path.join(db_name, f"{file_name}.lz4")


def _get_local_path(file_name, db_dir):
    """ Append local MIDAS_DB path to file(s) """
    if isinstance(file_name, list):
        return [os.path.join(db_dir, f"{fn}") for fn in file_name]
    return os.path.join(db_dir, f"{file_name}")


def _fetch_file_from_s3(packed_args):
    """ Fetch local path from AWS S3 dest path """
    s3_path, local_path = packed_args
    local_dir = os.path.dirname(local_path)
    if not os.path.isdir(local_dir):
        command(f"mkdir -p {local_dir}")
    if os.path.exists(local_path):
        return local_path
    return download_reference(s3_path, local_dir)
