#!/usr/bin/env python3
import os
import sys

from midas2.common.argparser import add_subcommand
from midas2.common.utils import tsprint, command, multithreading_map, num_physical_cores, pythonpath, InputStream
from midas2.common.utilities import decode_species_arg
from midas2.models.midasdb import MIDAS_DB
from midas2.params.inputs import MIDASDB_NAMES, MIDASDB_STATS


def list_midasdb(args):
    for dbname in MIDASDB_STATS.keys():
        dbdict = MIDASDB_STATS[dbname]
        nspecies = dbdict["species"]
        ngenomes = dbdict["genomes"]
        dbversion = dbdict["version"]
        print(f"{dbname} {ngenomes} genomes from {nspecies} species {dbversion}")


def init_midasdb(args):
    midasdb = MIDAS_DB(args.midasdb_dir, args.midasdb_name, 1)
    tsprint(f"Initate MIDASDB {args.midasdb_name} to {args.midasdb_dir}::start")
    midasdb.fetch_files("markerdb")
    midasdb.fetch_files("markerdb_models")
    midasdb.fetch_files("metadata")
    tsprint(f"Initate MIDASDB {args.midasdb_name} to {args.midasdb_dir}::finish")
    tsprint(f"Download chunks for MIDASDB {args.midasdb_name} to {args.midasdb_dir}::start")
    midasdb.fetch_files("chunks")
    tsprint(f"Download chunks for MIDASDB {args.midasdb_name} to {args.midasdb_dir}::finish")
    return True


def download_midasdb(args):
    assert args.species is not None or args.species_list is not None, f"Need to provide --species or --species_list for download task."
    if args.zzz_worker_mode:
        download_midasdb_worker(args)
    else:
        download_midasdb_master(args)
    return True


def download_midasdb_master(args):
    def sliced_work(i):
        n = args.num_cores
        worker_cmd = f"PYTHONPATH={pythonpath()} {sys.executable} -m midas2 database --download -s {i}:{n} --midasdb_name {args.midasdb_name} --midasdb_dir {args.midasdb_dir} --zzz_worker_mode {'--debug' if args.debug else ''}"
        try:
            command(worker_cmd)
        except Exception as error:
            if not args.debug:
                tsprint("Deleting untrustworthy dowdloaded databases. Specify --debug flag to keep.")
                command(f"rm -r {args.midasdb_dir}")
            raise error

    if  args.species == "all":
        multithreading_map(sliced_work, range(0, args.num_cores), num_threads=args.num_cores)
        return True

    midasdb = MIDAS_DB(os.path.abspath(args.midasdb_dir), args.midasdb_name, args.num_cores)
    species = midasdb.uhgg.species
    if args.species_list is not None:
        assert os.path.exists(args.species_list), f"Need to provide valid --species_list file."
        species_id_list = []
        with InputStream(args.species_list) as stream:
            for line in stream:
                species_id_list.append(line.strip())
    else:
        species_id_list = decode_species_arg(args, species)

    nspecies = len(species_id_list)
    tsprint(f"  Downloading MIDAS database for {nspecies} species::start")
    midasdb.fetch_files("repgenome", species_id_list)
    midasdb.fetch_files("pangenome", species_id_list)
    tsprint(f"  Downloading MIDAS database for {nspecies} species::finish")
    return True


def download_midasdb_worker(args):

    violation = "Please do not call download_midasdb_worker directly.  Violation"
    assert args.zzz_worker_mode, f"{violation}:  Missing --zzz_worker_mode arg."

    midasdb = MIDAS_DB(os.path.abspath(args.midasdb_dir), args.midasdb_name, 4) #<--- 4 way concurrency
    species = midasdb.uhgg.representatives

    species_id_list = decode_species_arg(args, species)
    i, n = args.species.split(":")

    tsprint(f"  Downloading MIDAS database for sliced species {i} with {n} cores in total::start")
    midasdb.fetch_files("repgenome", species_id_list)
    midasdb.fetch_files("pangenome", species_id_list)
    tsprint(f"  Downloading MIDAS database for sliced species {i} with {n} cores in total::finish")
    return True


def register_args(main_func):
    subparser = add_subcommand('database', main_func, help='List and download MIDAS DB from S3 to local')
    subparser.add_argument('--list',
                           action='store_true',
                           default=False,
                           help=f"List available MIDAS databases on S3.")
    subparser.add_argument('--init',
                           action='store_true',
                           default=False,
                           help=f"Initiate the download of MIDASDB with minimal files needed.")
    subparser.add_argument('--download',
                           action='store_true',
                           default=False,
                           help=f"Download select MIDAS database from S3 to local file.")
    subparser.add_argument('--midasdb_name',
                           dest='midasdb_name',
                           type=str,
                           default="uhgg",
                           choices=MIDASDB_NAMES,
                           help=f"MIDAS Database name.")
    subparser.add_argument('--midasdb_dir',
                           dest='midasdb_dir',
                           type=str,
                           default=".",
                           help=f"Path to local MIDAS Database.")
    subparser.add_argument('--num_cores',
                           dest='num_cores',
                           type=int,
                           metavar="INT",
                           default=num_physical_cores,
                           help=f"Number of physical cores to use ({num_physical_cores})")
    subparser.add_argument('-s',
                           '--species',
                           dest='species',
                           required=False,
                           help="species[,species...] whose MIDASDB(s) to build;  alternatively, species slice in format idx:modulus, e.g. 1:30, meaning build species whose ids are 1 mod 30; or, the special keyword 'all' meaning all species")
    subparser.add_argument('--species_list',
                           dest='species_list',
                           type=str,
                           required=False,
                           help="Path to list of species TXT file.")
    return main_func


@register_args
def main(args):
    #tsprint(f"Executing midas2 subcommand {args.subcommand} with args {vars(args)}.")
    if args.list:
        list_midasdb(args)
    if args.init:
        init_midasdb(args)
    if args.download:
        download_midasdb(args)
