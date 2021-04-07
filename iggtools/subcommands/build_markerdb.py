#!/usr/bin/env python3
import os
from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, retry, command, multithreading_map, find_files, upload, num_physical_cores, split, upload_star
from iggtools.models.midasdb import MIDAS_DB, MARKER_FILE_EXTS


@retry
def find_files_with_retry(f):
    return find_files(f)


def build_markerdb(args):
    """ Collate marker genes of repgenomes into phyeco.fa and phyeco.map """
    midas_db = MIDAS_DB(args.midas_db if args.midas_db else ".", num_cores=num_physical_cores)
    species = midas_db.uhgg.species
    do_upload = not args.midas_db

    build_marker_log_s3 = midas_db.get_target_layout("build_marker_log", remote=True)
    msg = f"Collating marker genes sequences."
    if find_files_with_retry(build_marker_log_s3):
        if not args.force and do_upload:
            tsprint(f"Destination {build_marker_log_s3} already exists.  Specify --force to overwrite.")
            return
        msg = msg.replace(msg.split(" ")[0], "Re-" + msg.split(" ")[0])
    tsprint(msg)

    build_marker_log = midas_db.get_target_layout("build_marker_log", remote=False)
    build_marker_subdir = os.path.dirname(build_marker_log)

    if not os.path.isdir(build_marker_subdir):
        command(f"mkdir -p {build_marker_subdir}")

    with open(build_marker_log, "w") as slog:
        slog.write(msg + "\n")

    # Fetch marker genes fasta-file and map-file from s3
    list_of_repgenome_marker_seqs = midas_db.fetch_files("marker_genes_fa", species.keys())
    list_of_repgenome_marker_maps = midas_db.fetch_files("marker_genes_map", species.keys())

    # Collate to phyeco.fa and phyeco.map
    phyeco_seqs = midas_db.get_target_layout("marker_db", remote=False, component="fa")
    for _ in split(list_of_repgenome_marker_seqs.values(), 20):
        command("cat " + " ".join(_) + f" >> {phyeco_seqs}")

    phyeco_maps = midas_db.get_target_layout("marker_db", remote=False, component="map")
    for _ in split(list_of_repgenome_marker_maps.values(), 20):
        command("cat " + " ".join(_) + f" >> {phyeco_maps}")

    # Build hs-blastn index for the collated phyeco sequences
    cmd_index = f"hs-blastn index {phyeco_seqs} &>> {build_marker_log}"
    with open(f"{build_marker_log}", "a") as slog:
        slog.write(cmd_index + "\n")
    command(cmd_index)

    # Upload generated fasta and index files
    if do_upload:
        upload_tasks = []
        for ext in MARKER_FILE_EXTS:
            local_file = midas_db.get_target_layout("marker_db", remote=False, component=ext)
            s3_file = midas_db.get_target_layout("marker_db", remote=True, component=ext)
            upload_tasks.append((local_file, s3_file))
        multithreading_map(upload_star, upload_tasks)
        upload(build_marker_log, build_marker_log_s3, check=False) # Upload LOG file last as a checkout point


def register_args(main_func):
    subparser = add_subcommand('build_markerdb', main_func, help='collate marker genes for repgresentative genomes')
    subparser.add_argument('--midas_db', # when provided, user build customer phyeco database only locally
                           dest='midas_db',
                           type=str,
                           required=False,
                           metavar="CHAR",
                           help=f"local MIDAS DB which mirrors the s3 IGG db")
    return main_func


@register_args
def main(args):
    tsprint(f"Executing iggtools subcommand {args.subcommand} with args {vars(args)}.")
    build_markerdb(args)
