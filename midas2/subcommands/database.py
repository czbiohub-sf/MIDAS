#!/usr/bin/env python3
import os
import sys

from midas2.common.argparser import add_subcommand
from midas2.common.utils import tsprint, command, multithreading_map, num_physical_cores, pythonpath
from midas2.common.utilities import decode_species_arg
from midas2.models.midasdb import MIDAS_DB
from midas2.params.inputs import MIDASDB_NAMES, MIDASDB_VERSION


def list_midasdb(args):
    for dbname in MIDASDB_NAMES:
        if dbname == "testdb":
            continue
        midasdb = MIDAS_DB(args.midasdb_dir, dbname, 1)
        nspecies = len(midasdb.uhgg.species)
        ngenomes = len(midasdb.uhgg.genomes)
        dbversion = MIDASDB_VERSION[dbname]
        print(f"{dbname} {ngenomes} genomes from {nspecies} species {dbversion}")


def download_midasdb(args):
    midasdb = MIDAS_DB(args.midasdb_dir, args.midasdb_name, 1)
    midasdb.fetch_files("marker_db")
    midasdb.fetch_files("marker_db_hmm_cutoffs")

    if args.markerdb_only:
        return True

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

    multithreading_map(sliced_work, range(0, args.num_cores), num_threads=args.num_cores)


def download_midasdb_worker(args):

    violation = "Please do not call download_midasdb_worker directly.  Violation"
    assert args.zzz_worker_mode, f"{violation}:  Missing --zzz_worker_mode arg."

    midasdb = MIDAS_DB(os.path.abspath(args.midasdb_dir), args.midasdb_name, 4) #<--- 4 way concurrency
    species = midasdb.uhgg.representatives

    species_id_list = decode_species_arg(args, species)
    i, n = args.species.split(":")
    snps_chunk_size = 1000000
    genes_chunk_size = 50000

    tsprint(f"  Downloading MIDAS database for sliced species {i} with {n} cores in total::start")
    midasdb.fetch_files("annotation_fna", species_id_list)
    midasdb.fetch_files("annotation_genes", species_id_list)
    midasdb.fetch_files("annotation_ffn", species_id_list)
    midasdb.fetch_files("pangenome_centroids", species_id_list)
    midasdb.fetch_files("pangenome_cluster_info", species_id_list)
    midasdb.fetch_chunks("chunks_sites_run", species_id_list, snps_chunk_size)
    midasdb.fetch_chunks("chunks_sites_merge", species_id_list, snps_chunk_size)
    midasdb.fetch_chunks("chunks_centroids", species_id_list, genes_chunk_size)
    tsprint(f"  Downloading MIDAS database for sliced species {i} with {n} cores in total::finish")
    return True


def register_args(main_func):
    subparser = add_subcommand('database', main_func, help='List and download MIDAS DB from S3 to local')
    subparser.add_argument('--list',
                           action='store_true',
                           default=False,
                           help=f"List available MIDAS databases on S3.")
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
                           help=f"Local MIDAS Database path mirroing S3.")
    subparser.add_argument('--markerdb_only',
                           action='store_true',
                           default=False,
                           help=f"Only download marker databases related files.")
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
                           help="species slice in format idx:modulus, e.g. 1:30.")
    return main_func


@register_args
def main(args):
    #tsprint(f"Executing midas2 subcommand {args.subcommand} with args {vars(args)}.")
    if args.list:
        list_midasdb(args)
    if args.download:
        download_midasdb(args)
