import os
from collections import defaultdict
from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, InputStream, OutputStream, select_from_tsv, retry, command, multithreading_map, find_files, upload, num_physical_cores, split, upload_star, download_reference
from iggtools.models.uhgg import MIDAS_IGGDB, MARKER_FILE_EXTS, get_uhgg_layout
from iggtools.params import inputs, outputs
from iggtools.params.schemas import MARKER_INFO_SCHEMA, PAN_GENE_INFO_SCHEMA


@retry
def find_files_with_retry(f):
    return find_files(f)


# Find out which pan-gene is the marker-gene. This is a UHGG related add on scripts.
def map_marker_to_centroids(args):
    """ Identify which gene clusters have repgenome's marker genes """
    midas_iggdb = MIDAS_IGGDB(args.midas_iggdb)
    representatives = midas_iggdb.uhgg.representatives

    # Alternatively, you can also fetch the collated phyeco.map for all species
    fetched_marker_genes_mapfile = midas_iggdb.fetch_files(get_uhgg_layout(species_id="", component="map")["marker_db"])

    log_remote = midas_iggdb.get_target_layout("marker_centroids_log", remote=True)
    msg = f"Finding centroids genes for marker genes."
    if find_files_with_retry(log_remote):
        if not args.force:
            tsprint(f"Destination {log_remote} already exists.  Specify --force to overwrite.")
            return
        msg = msg.replace(msg.split(" ")[0], "Re-" + msg.split(" ")[0])
    tsprint(msg)

    log_local = midas_iggdb.get_target_layout("marker_centroids_log", remote=False)
    subdir = os.path.dirname(log_local)
    if not os.path.isdir(subdir):
        command(f"mkdir -p {subdir}")
    with open(log_local, "w") as slog:
        slog.write(msg + "\n")

    def species_work(species_id):
        # Fetch per-species marker_genes_mapfile
        fetched_marker_map = midas_iggdb.fetch_files(get_uhgg_layout(species_id, "markers.map", representatives[species_id])["marker_genes"])
        fetched_gene_info = midas_iggdb.fetch_files(get_uhgg_layout(species_id, "gene_info.txt")["pangenome_file"])

        # Parse the marker gene id of the representative genomes from the phyeco.mapfile. Refer to collate_repgenome_markers.
        markers = dict()
        awk_command = f"awk \'$1 == \"{species_id}\"\'"
        with InputStream(fetched_marker_map, awk_command) as stream:
            for gene_id, marker_id in select_from_tsv(stream, ["gene_id", "marker_id"], schema=MARKER_INFO_SCHEMA):
                assert marker_id not in markers, f"marker {marker_id} for species {species_id} corresponds to multiple gene_ids."
                markers[gene_id] = marker_id

        mc_file_local = midas_iggdb.get_target_layout("marker_centroids", remote=False, component="", species_id=species_id)
        mc_file_remote = midas_iggdb.get_target_layout("marker_centroids", remote=True, component="", species_id=species_id)

        marker_centroids_subdir = os.path.dirname(mc_file_local)
        if not os.path.isdir(os.path.dirname(mc_file_local)):
            command(f"mkdir -p {marker_centroids_subdir}")

        # Filter gene_info.txt by marker_ids
        with OutputStream(mc_file_local) as ostream:
            ostream.write("\t".join(["marker_id"] + list(PAN_GENE_INFO_SCHEMA.keys())) + "\n")
            with InputStream(fetched_gene_info) as stream:
                for row in select_from_tsv(stream, selected_columns=PAN_GENE_INFO_SCHEMA, result_structure=dict):
                    if row["gene_id"] in markers.keys():
                        ostream.write("\t".join([markers[row["gene_id"]]] + list(row.values())) + "\n")

        # Update to s3
        upload(mc_file_local, mc_file_remote, check=False)

        with open(f"{log_local}", "a") as slog:
            slog.write(f"Species {species_id} finished" + "\n")

    multithreading_map(species_work, midas_iggdb.uhgg.species.keys(), num_physical_cores)

    # Upload log_local log file in the last
    upload(collate_log, log_remote, check=False)


def collate_repgenome_markers(args):
    """ Collate marker genes of repgenomes into phyeco.fa and phyeco.map """
    midas_iggdb = MIDAS_IGGDB(args.midas_iggdb)
    species = midas_iggdb.uhgg.species

    collate_log_remote = midas_iggdb.get_target_layout("marker_collate_log", remote=True)
    msg = f"Collating marker genes sequences."
    if find_files_with_retry(collate_log_remote):
        if not args.force:
            tsprint(f"Destination {collate_log_remote} already exists.  Specify --force to overwrite.")
            return
        msg = msg.replace(msg.split(" ")[0], "Re-" + msg.split(" ")[0])
    tsprint(msg)

    collate_log = midas_iggdb.get_target_layout("marker_collate_log", remote=False)
    collate_subdir = os.path.dirname(collate_log)
    if not args.debug:
        command(f"rm -rf {collate_subdir}")
    if not os.path.isdir(collate_subdir):
        command(f"mkdir -p {collate_subdir}")
    with open(collate_log, "w") as slog:
        slog.write(msg + "\n")

    # Fetch marker genes fasta-file and map-file from s3
    marker_genes_fasta = midas_iggdb.fetch_files("marker_genes_fa", species.keys())
    marker_genes_maps = midas_iggdb.fetch_files("marker_genes_map", species.keys())

    # Collate to phyeco.fa and phyeco.map
    collated_genes_fa = midas_iggdb.get_target_layout("marker_db", remote=False, component="fa")
    for marker_fa_files in split(marker_genes_fasta.values(), 20):
        command("cat " + " ".join(marker_fa_files) + f" >> {collated_genes_fa}")

    collaged_genes_map = midas_iggdb.get_target_layout("marker_db", remote=False, component="map")
    for marker_map_files in split(marker_genes_maps.values(), 20):
        command("cat " + " ".join(marker_map_files) + f" >> {collaged_genes_map}")

    # Build hs-blastn index for the collated phyeco sequences
    #? do we have to provide the collage_subdir to hs-blastn?
    # can only check this one out on S3, wait until the marker_centroids works
    cmd_index = f"cd {collate_subdir}; hs-blastn index {collated_genes_fa} &>> {collate_log}"
    with open(f"{collate_log}", "a") as slog:
        slog.write(cmd_index + "\n")
    command(cmd_index)

    # Upload generated fasta and index files
    upload_tasks = []
    for ext in MARKER_FILE_EXTS:
        local_file = midas_iggdb.get_target_layout("marker_db", remote=False, component=ext)
        s3_file = midas_iggdb.get_target_layout("marker_db", remote=True, component=ext)
        upload_tasks.append((local_file, s3_file))
    multithreading_map(upload_star, upload_tasks)

    # Upload the log file in the last
    upload(collate_log, collate_log_remote, check=False)

    # Clean up
    if not args.debug:
        command(f"rm -rf {collate_subdir}", check=False)


def register_args(main_func):
    subparser = add_subcommand('collate_repgenome_markers', main_func, help='collate marker genes for repgresentative genomes')
    subparser.add_argument('--midas_iggdb',
                           dest='midas_iggdb',
                           type=str,
                           required=True,
                           metavar="CHAR",
                           help=f"local MIDAS DB which mirrors the s3 IGG db")
    subparser.add_argument('--collate_repgenome_markers',
                           action='store_true',
                           default=False,
                           help=f"collate marker genes and mapfile for the repgenomes")
    subparser.add_argument('--map_marker_to_centroids',
                           action='store_true',
                           default=False,
                           help=f"Identify the centroids genes for rep marker genes")
    return main_func


@register_args
def main(args):
    tsprint(f"Executing iggtools subcommand {args.subcommand} with args {vars(args)}.")
    if args.collate_repgenome_markers:
        collate_repgenome_markers(args)
    if args.map_marker_to_centroids:
        map_marker_to_centroids(args)
