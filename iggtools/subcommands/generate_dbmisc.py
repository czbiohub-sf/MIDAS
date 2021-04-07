#!/usr/bin/env python3
import os
import sys
from multiprocessing import Semaphore
import gffutils
from iggtools.common.argparser import add_subcommand, SUPPRESS
from iggtools.common.utils import tsprint, InputStream, OutputStream, download_reference, select_from_tsv, retry, command, multithreading_map, multiprocessing_map, find_files, upload, num_physical_cores, pythonpath
from iggtools.models.uhgg import UHGG, get_uhgg_layout, destpath
from iggtools.params import outputs
from iggtools.subcommands.import_uhgg import decode_genomes_arg
from iggtools.subcommands.build_pangenome import decode_species_arg
from iggtools.params.schemas import PAN_GENE_INFO_SCHEMA, MARKER_INFO_SCHEMA, PAN_GENE_LENGTH_SCHEMA, CLUSTER_INFO_SCHEMA


# Up to this many concurrent species builds.
CONCURRENT_BUILDS = Semaphore(num_physical_cores-1)


@retry
def find_files_with_retry(f):
    return find_files(f)


def convert_gff_to_genes(gff3_file, genes_file):
    """ Convert GFF3 features format into genes.feature """
    db = gffutils.create_db(gff3_file, f"{gff3_file}.db")
    with OutputStream(genes_file) as stream:
        stream.write("\t".join(["gene_id", "contig_id", "start", "end", "strand", "gene_type"]) + "\n")
        for feature in db.all_features():
            if feature.source == "prokka":
                continue
            if "ID" not in feature.attributes: #CRISPR
                continue
            seqid = feature.seqid
            start = feature.start
            stop = feature.stop
            strand = feature.strand
            gene_id = feature.attributes['ID'][0]
            locus_tag = feature.attributes['locus_tag'][0]
            assert gene_id == locus_tag
            gene_type = feature.featuretype
            stream.write("\t".join([gene_id, seqid, str(start), str(stop), strand, gene_type]) + "\n")
    return True


def generate_gene_features(args):
    if args.zzz_slave_toc:
        generate_gene_features_slave(args)
    else:
        generate_gene_features_master(args)


def generate_gene_features_master(args):

    local_toc = download_reference(outputs.genomes)

    db = UHGG(local_toc)
    species_for_genome = db.genomes

    def genome_work(genome_id):
        assert genome_id in species_for_genome, f"Genome {genome_id} is not in the database."
        species_id = species_for_genome[genome_id]

        dest_file = destpath(get_uhgg_layout(species_id, "genes", genome_id)["annotation_file"])
        msg = f"Builing gene features for genome {genome_id} from species {species_id}."
        if find_files_with_retry(dest_file):
            if not args.force:
                tsprint(f"Destination {dest_file} for genome {genome_id} gene features already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Importing", "Reimporting")

        with CONCURRENT_BUILDS:
            tsprint(msg)
            logfile = get_uhgg_layout(species_id, "", genome_id)["gene_features_log"]
            slave_log = os.path.basename(logfile)
            slave_subdir = f"misc/{species_id}__{genome_id}"
            if not args.debug:
                command(f"rm -rf {slave_subdir}")
            if not os.path.isdir(slave_subdir):
                command(f"mkdir -p {slave_subdir}")

            # Recurisve call via subcommand.  Use subdir, redirect logs.
            slave_cmd = f"cd {slave_subdir}; PYTHONPATH={pythonpath()} {sys.executable} -m iggtools generate_dbmisc --generate_gene_features --genome {genome_id} --zzz_slave_mode --zzz_slave_toc {os.path.abspath(local_toc)} {'--debug' if args.debug else ''} &>> {slave_log}"
            with open(f"{slave_subdir}/{slave_log}", "w") as slog:
                slog.write(msg + "\n")
                slog.write(slave_cmd + "\n")

            try:
                command(slave_cmd)
            finally:
                upload(f"{slave_subdir}/{slave_log}", destpath(logfile), check=False)
                if not args.debug:
                    command(f"rm -rf {slave_subdir}", check=False)

    genome_id_list = decode_genomes_arg(args, species_for_genome)
    multithreading_map(genome_work, genome_id_list, num_physical_cores)


def generate_gene_features_slave(args):
    """
    https://github.com/czbiohub/iggtools/wiki
    """

    violation = "Please do not call generate_gene_features_slave directly.  Violation"
    assert args.zzz_slave_mode, f"{violation}:  Missing --zzz_slave_mode arg."
    assert os.path.isfile(args.zzz_slave_toc), f"{violation}: File does not exist: {args.zzz_slave_toc}"

    db = UHGG(args.zzz_slave_toc)
    species_for_genome = db.genomes

    genome_id = args.genomes
    species_id = species_for_genome[genome_id]

    input_annot_file = f"{genome_id}.gff"
    download_reference(destpath(get_uhgg_layout(species_id, "gff", genome_id)["annotation_file"]))

    last_output = f"{genome_id}.genes"
    dest_file = destpath(get_uhgg_layout(species_id, "genes", genome_id)["annotation_file"])
    assert convert_gff_to_genes(input_annot_file, last_output)

    upload(last_output, dest_file)


def parse_markers_per_genome(packed_ids):
    """ Extract makrer_id, gene_id pairs from the MAP file per genome """
    species_id, genome_id = packed_ids
    marker_by_genome = dict()
    with InputStream(destpath(get_uhgg_layout(species_id, "markers.map", genome_id)["marker_genes"])) as stream:
        for gene_id, marker_id in select_from_tsv(stream, ["gene_id", "marker_id"], schema=MARKER_INFO_SCHEMA):
            marker_by_genome[gene_id] = marker_id
    return marker_by_genome


def join_marker_centroids(args):
    if args.zzz_slave_toc:
        join_marker_centroids_slave(args)
    else:
        join_marker_centroids_master(args)


def join_marker_centroids_master(args):

    local_toc = download_reference(outputs.genomes)

    db = UHGG(local_toc)
    species = db.species

    def species_work(species_id):
        assert species_id in species, f"Species {species_id} is not in the database."
        species_genomes = species[species_id]

        dest_file = destpath(get_uhgg_layout(species_id, "txt")["marker_centroids"])
        msg = f"Building marker genes to centroids for species {species_id} with {len(species_genomes)} total genomes."
        if find_files_with_retry(dest_file):
            if not args.force:
                tsprint(f"Destination {dest_file} for species {species_id} pangenome already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Building", "Rebuilding")

        with CONCURRENT_BUILDS:
            tsprint(msg)
            logfile = get_uhgg_layout(species_id)["cluster_info_log"]
            slave_log = os.path.basename(logfile)
            slave_subdir = f"marker_centroids/{species_id}"
            if not args.debug:
                command(f"rm -rf {slave_subdir}")
            if not os.path.isdir(slave_subdir):
                command(f"mkdir -p {slave_subdir}")

            # Recurisve call via subcommand.  Use subdir, redirect logs.
            slave_cmd = f"cd {slave_subdir}; PYTHONPATH={pythonpath()} {sys.executable} -m iggtools generate_dbmisc --join_marker_centroids -s {species_id} --zzz_slave_mode --zzz_slave_toc {os.path.abspath(local_toc)} {'--debug' if args.debug else ''} &>> {slave_log}"
            with open(f"{slave_subdir}/{slave_log}", "w") as slog:
                slog.write(msg + "\n")
                slog.write(slave_cmd + "\n")

            try:
                command(slave_cmd)
            finally:
                upload(f"{slave_subdir}/{slave_log}", destpath(get_uhgg_layout(species_id, "log")["marker_centroids"]), check=False)
                if not args.debug:
                    command(f"rm -rf {slave_subdir}", check=False)

    species_id_list = decode_species_arg(args, species)
    multithreading_map(species_work, species_id_list, num_physical_cores)


def join_marker_centroids_slave(args):
    violation = "Please do not call join_marker_centroids_slave directly.  Violation"
    assert args.zzz_slave_mode, f"{violation}:  Missing --zzz_slave_mode arg."
    assert os.path.isfile(args.zzz_slave_toc), f"{violation}: File does not exist: {args.zzz_slave_toc}"

    db = UHGG(args.zzz_slave_toc)
    species = db.species
    species_id = args.species

    assert species_id in species, f"{violation}: Species {species_id} is not in the database."

    species_genomes = species[species_id]
    species_genomes_ids = species_genomes.keys()

    dest_file = destpath(get_uhgg_layout(species_id, "txt")["marker_centroids"])
    local_file = f"{species_id}.txt"

    # Read in all genome's phyeco.mapfile
    list_marker_dicts = multiprocessing_map(parse_markers_per_genome, ((species_id, genome_id) for genome_id in species_genomes_ids))
    # flatten the list of dict
    dict_of_markers = {k: v for md in list_marker_dicts for k, v, in md.items()}

    # Read in species/genomes/mapfile
    gene_info_file = destpath(get_uhgg_layout(species_id, "gene_info.txt")["pangenome_file"])

    # Filter gene_info.txt by marker_ids and write to file
    with OutputStream(local_file) as ostream:
        ostream.write("\t".join(["marker_id"] + list(PAN_GENE_INFO_SCHEMA.keys())) + "\n")
        with InputStream(gene_info_file) as stream:
            # Loop over the ALl the genes, locate those that are marker genes, and write the
            for row in select_from_tsv(stream, selected_columns=PAN_GENE_INFO_SCHEMA, result_structure=dict):
                if row["gene_id"] in dict_of_markers.keys():
                    ostream.write("\t".join([dict_of_markers[row["gene_id"]]] + list(row.values())) + "\n")

    # Update to s3
    command(f"aws s3 rm --recursive {dest_file}")
    upload(local_file, dest_file, check=False)




def species_cluster_work(args):

    species_id, args_force, args_debug = args
    dest_file = destpath(get_uhgg_layout(species_id, "cluster_info.txt")["pangenome_file"])

    msg = f"Building cluster info file for species {species_id}."
    if find_files_with_retry(dest_file):
        if not args_force:
            tsprint(f"Destination {dest_file} for species {species_id} pangenome already exists.  Specify --force to overwrite.")
            return
        msg = msg.replace("Building", "Rebuilding")

    tsprint(msg)
    logfile = get_uhgg_layout(species_id)["cluster_info_log"]
    slave_log = os.path.basename(logfile)
    slave_subdir = f"cluster_info/{species_id}"
    if not args_force:
        command(f"rm -rf {slave_subdir}")
    if not os.path.isdir(slave_subdir):
        command(f"mkdir -p {slave_subdir}")


    gene_info_file = destpath(get_uhgg_layout(species_id, "gene_info.txt")["pangenome_file"])
    centroids_dict = {}
    cols = list(PAN_GENE_INFO_SCHEMA.keys())[1:]
    with InputStream(gene_info_file) as stream:
        for r in select_from_tsv(stream, selected_columns=cols, result_structure=dict):
            if r["centroid_99"] not in centroids_dict:
                centroids_dict[r["centroid_99"]] = r

    gene_length_file = destpath(get_uhgg_layout(species_id, "genes.len")["pangenome_file"])
    gene_length_dict = dict()
    with InputStream(gene_length_file) as stream:
        for r in select_from_tsv(stream, schema = PAN_GENE_LENGTH_SCHEMA, result_structure=dict):
            gene_length_dict[r["gene_id"]] = r["gene_length"]


    last_output = f"{slave_subdir}/cluster_info.txt"
    with OutputStream(last_output) as stream:
        stream.write("\t".join(CLUSTER_INFO_SCHEMA.keys()) + "\n")
        for record in centroids_dict.values():
            gene_len = gene_length_dict[record["centroid_99"]]
            val = list(record.values())
            val.append(gene_len)
            stream.write("\t".join(map(str, val)) + "\n")

    with open(f"{slave_subdir}/{slave_log}", "w") as slog:
        slog.write(msg + "\n")

    upload(last_output, dest_file)
    upload(f"{slave_subdir}/{slave_log}", destpath(logfile), check=False)
    if not args_debug:
        command(f"rm -rf {slave_subdir}", check=False)


def generate_cluster_info(args):
    local_toc = download_reference(outputs.genomes)

    db = UHGG(local_toc)
    species = db.species

    species_id_list = decode_species_arg(args, species)
    args_list = [(species_id, args.force, args.debug) for species_id in species_id_list]
    multiprocessing_map(species_cluster_work, args_list, num_physical_cores)


def register_args(main_func):
    subparser = add_subcommand('generate_dbmisc', main_func, help='Generate variety of db related files desired by MIDAS')
    subparser.add_argument('--genomes',
                           dest='genomes',
                           required=False,
                           help="genome[,genome...] to import;  alternatively, slice in format idx:modulus, e.g. 1:30, meaning annotate genomes whose ids are 1 mod 30; or, the special keyword 'all' meaning all genomes")
    subparser.add_argument('-s',
                           '--species',
                           dest='species',
                           required=False,
                           help="species[,species...] whose pangenome(s) to build;  alternatively, species slice in format idx:modulus, e.g. 1:30, meaning build species whose ids are 1 mod 30; or, the special keyword 'all' meaning all species")
    subparser.add_argument('--generate_gene_features',
                           action='store_true',
                           default=False,
                           help=f"Generate gene features for each genomes")
    subparser.add_argument('--join_marker_centroids',
                           action='store_true',
                           default=False,
                           help=f"Join marker genes with centroid_99 genes for each species")
    subparser.add_argument('--generate_cluster_info',
                           action='store_true',
                           default=False,
                           help=f"Generate cluster_info.txt used in midas_merge_genes.")
    subparser.add_argument('--zzz_slave_toc',
                           dest='zzz_slave_toc',
                           required=False,
                           help=SUPPRESS) # "reserved to pass table of contents from master to slave"
    return main_func


@register_args
def main(args):
    tsprint(f"Executing iggtools subcommand {args.subcommand} with args {vars(args)}.")
    if args.generate_gene_features:
        generate_gene_features(args)
    if args.join_marker_centroids:
        join_marker_centroids(args)
    if args.generate_cluster_info:
        generate_cluster_info(args)
