#!/usr/bin/env python3
from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, command
from iggtools.models.midasdb import MIDAS_DB
from iggtools.models.species import parse_species

CONCURRENT_DOWNLOAD = 20

def register_args(main_func):
    subparser = add_subcommand('download_midasdb', main_func, help='download midas database from S3 to local')
    subparser.add_argument('--midasdb_name',
                           dest='midasdb_name',
                           type=str,
                           default="uhgg",
                           choices=['uhgg', 'gtdb', 'testdb'],
                           help=f"MIDAS Database name.")
    subparser.add_argument('--midasdb_dir',
                           dest='midasdb_dir',
                           type=str,
                           default=".",
                           help=f"Local MIDAS Database path mirroing S3.")
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
                           default=CONCURRENT_DOWNLOAD,
                           help=f"Number of physical cores to use ({CONCURRENT_DOWNLOAD})")
    subparser.add_argument('--species_list',
                           dest='species_list',
                           type=str,
                           metavar="CHAR",
                           help=f"Path to file OR comma separated list of species ids")
    return main_func


def download_midasdb(args):
    try:
        midas_db = MIDAS_DB(args.midasdb_dir, args.midasdb_name, args.num_cores)
        species = midas_db.uhgg.representatives
        species_all = set(species.keys())

        species_list = set(parse_species(args))

        if args.db_type == "basic":
            list_of_species = list(species_list & species_all)
            list_of_species = list(species_all)[:5] if not list_of_species else list_of_species
        else:
            list_of_species = list(species_all)
        tsprint(len(list_of_species))

        midas_db.fetch_files("representative_genome", list_of_species)
        midas_db.fetch_files("pangenome_centroids", list_of_species)
        midas_db.fetch_files("pangenome_cluster_info", list_of_species)

        midas_db.fetch_files("marker_db")
        midas_db.fetch_files("marker_db_hmm_cutoffs")

    except Exception as error:
        if not args.debug:
            tsprint("Deleting untrustworthy dowdloaded databases. Specify --debug flag to keep.")
            command(f"rm -r {args.midasdb_dir}")
        raise error


@register_args
def main(args):
    tsprint(f"Executing iggtools subcommand {args.subcommand} with args {vars(args)}.")
    download_midasdb(args)
