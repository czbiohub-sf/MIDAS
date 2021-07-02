#!/usr/bin/env python3
import os
from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, retry, InputStream, OutputStream, command, multithreading_map, find_files, upload, num_physical_cores, split, upload_star
from iggtools.models.midasdb import MIDAS_DB, MARKER_FILE_EXTS

CONCURRENT_DOWNLOAD = 20

def register_args(main_func):
    subparser = add_subcommand('download_midasdb', main_func, help='download midas database from S3 to local')
    subparser.add_argument('--midas_db', # when provided, user build customer phyeco database only locally
                           dest='midas_db',
                           type=str,
                           required=True,
                           metavar="CHAR",
                           help=f"local MIDAS DB which mirrors the s3 IGG db")
    subparser.add_argument('--db_type', # when provided, user build customer phyeco database only locally
                           dest='db_type',
                           type=str,
                           required=False,
                           metavar="CHAR",
                           default="basic",
                           help=f"Download basic version or full db (basic recommanded).")
    subparser.add_argument('--num_cores',
                           dest='num_cores',
                           type=int,
                           metavar="INT",
                           default=num_physical_cores,
                           help=f"Number of physical cores to use ({num_physical_cores})")
    subparser.add_argument('--species_list',
                           dest='species_list',
                           type=str,
                           metavar="CHAR",
                           help=f"Comma separated list of species ids")
    subparser.add_argument('--species_file',
                           dest='species_file',
                           type=str,
                           metavar="CHAR",
                           help=f"Path to list of species ids")
    return main_func


def scan_species(filename):
    splist = []
    with InputStream(filename) as stream:
        for line in stream:
            splist.append(line.strip("\n"))
    return splist


def download_midasdb(args):
    try:
        midas_db = MIDAS_DB(args.midas_db, args.num_cores)

        species_list1 = args.species_list.split(",") if args.species_list else []
        species_list2 = scan_species(args.species_file) if args.species_file else []
        species_list = set(species_list1 + species_list2)

        species = midas_db.uhgg.representatives
        species_all = set(species.keys())

        if args.db_type == "basic":
            list_of_species = list(species_list & species_all)
            list_of_species = list(species_all)[:5] if not list_of_species else list_of_species
        else:
            list_of_species = list(species_all)

        tsprint(list_of_species)

        # Marker genes related database files
        midas_db.fetch_files("marker_centroids", list_of_species)
        midas_db.fetch_files("cluster_info", list_of_species)

        midas_db.fetch_files("prokka_genome", list_of_species)
        midas_db.fetch_files("centroids", list_of_species)

        marker_db_files = midas_db.fetch_files("marker_db")
        marker_db_hmm_cutoffs = midas_db.fetch_files("marker_db_hmm_cutoffs")
        
    except Exception as error:
        if not args.debug:
            tsprint("Deleting untrustworthy dowdloaded databases. Specify --debug flag to keep.")
            command(f"rm -r {args.midas_db}")
        raise error


@register_args
def main(args):
    tsprint(f"Executing iggtools subcommand {args.subcommand} with args {vars(args)}.")
    download_midasdb(args)
