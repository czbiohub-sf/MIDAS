#!/usr/bin/env python3
import json
import os
from midas2.common.argparser import add_subcommand
from midas2.common.utils import tsprint, num_physical_cores, command
from midas2.common.bowtie2 import build_bowtie2_db
from midas2.models.midasdb import MIDAS_DB
from midas2.models.species import parse_species, filter_species
from midas2.params.inputs import MIDASDB_NAMES


DEFAULT_SAMPLE_COUNTS = 2
DEFAULT_PRUNE_CUTOFF = 0.4

def register_args(main_func):
    subparser = add_subcommand('build_bowtie2db', main_func, help='Build rep-genome and pan-genome Bowtie2 indexes for list of species')

    subparser.add_argument('--bt2_indexes_dir',
                           dest='bt2_indexes_dir',
                           type=str,
                           metavar="CHAR",
                           required=True,
                           help="Path to bowtie2 indexes directory")
    subparser.add_argument('--bt2_indexes_name',
                           dest='bt2_indexes_name',
                           type=str,
                           metavar="CHAR",
                           default="repgenomes",
                           choices=['repgenomes', 'pangenomes'],
                           help="Bowtie2 db name to build")

    subparser.add_argument('--midasdb_name',
                           dest='midasdb_name',
                           type=str,
                           default="uhgg",
                           choices=MIDASDB_NAMES,
                           help="MIDAS Database name.")
    subparser.add_argument('--midasdb_dir',
                           dest='midasdb_dir',
                           type=str,
                           default=".",
                           help="Path to local MIDAS Database.")

    subparser.add_argument('--species_list',
                           dest='species_list',
                           type=str,
                           metavar="CHAR",
                           help="Comma separated list of species ids OR path to species list txt file")
    subparser.add_argument('--species_profile',
                           dest='species_profile',
                           type=str,
                           metavar="CHAR",
                           help="Path to species coverage TSV file with species_id")
    subparser.add_argument('--select_by',
                           dest='select_by',
                           type=str,
                           metavar="CHAR",
                           default="sample_counts",
                           help="Column from species_prevalence based on which to select species (sample_counts).")
    subparser.add_argument('--select_threshold',
                           dest='select_threshold',
                           type=str,
                           metavar="CHAR",
                           default=str(DEFAULT_SAMPLE_COUNTS),
                           help="Comman separated correponsding cutoff to select_by (>XX) ({DEFAULT_SAMPLE_COUNTS}, )")

    subparser.add_argument('--num_cores',
                           dest='num_cores',
                           type=int,
                           metavar="INT",
                           default=num_physical_cores,
                           help=f"Number of physical cores to use ({num_physical_cores})")
    # Prune centroids_99 for cleaner reads mapping
    subparser.add_argument('--prune_centroids',
                           action='store_true',
                           default=False,
                           help='Prune shorter centroids99 within each centroids95 cluster')
    subparser.add_argument('--prune_method',
                           dest='prune_method',
                           type=str,
                           default='max',
                           choices=['max', 'median', 'mean'],
                           help="Prune methods: max, median or mean.")
    subparser.add_argument('--prune_cutoff',
                           dest='prune_cutoff',
                           type=float,
                           default=DEFAULT_PRUNE_CUTOFF,
                           help=f"Prune cutoff: for each centroid_95, centroid99 shorter than {DEFAULT_PRUNE_CUTOFF} of the chosen method are pruned for reading mapping.")
    subparser.add_argument('--remove_singleton',
                           action='store_true',
                           default=False,
                           help='Remove c75 clusters with only one gene member in species with more genomes.')
    return main_func


def fetch_pruned_targets(midas_db, species_of_interest, opts, cutoff, remove_singleton):
    if remove_singleton:
        return {spid: midas_db.get_target_layout("pruned_centroids_rs", False, spid, opts, cutoff) for spid in species_of_interest}
    return {spid: midas_db.get_target_layout("pruned_centroids", False, spid, opts, cutoff) for spid in species_of_interest}


def build_bowtie2db(args):

    try:
        if args.species_list:
            species_of_interest = parse_species(args)
        elif args.species_profile and args.select_by and args.select_threshold:
            species_of_interest = filter_species(args.species_profile, args.select_by, args.select_threshold)
        else:
            raise Exception("Need to provide either species_list or species_profile as input arguments")
        tsprint(f"MIDAS2::build_bowtie2db::build bt2 indexees for the listed species: {species_of_interest}")

        # Fetch MIDAS Reference Database Files
        midas_db = MIDAS_DB(os.path.abspath(args.midasdb_dir), args.midasdb_name, args.num_cores)

        if args.bt2_indexes_name == "repgenomes":
            tsprint("MIDAS2::build_bowtie2_repgenomes_indexes::start")
            midas_db.fetch_files("repgenome", species_of_interest)
            contigs_files = midas_db.fetch_files("representative_genome", species_of_interest)
            build_bowtie2_db(args.bt2_indexes_dir, args.bt2_indexes_name, contigs_files, args.num_cores)
            tsprint("MIDAS2::build_bowtie2_repgenomes_indexes::finish")

        if args.bt2_indexes_name == "pangenomes":
            tsprint("MIDAS2::build_bowtie2_pangenomes_indexes::start")
            midas_db.fetch_files("pangenome", species_of_interest)
            if args.prune_centroids:
                centroids_files = fetch_pruned_targets(midas_db, species_of_interest, args.prune_method, args.prune_cutoff, args.remove_singleton)
            else:
                centroids_files = midas_db.fetch_files("pangenome_centroids", species_of_interest)
            build_bowtie2_db(args.bt2_indexes_dir, args.bt2_indexes_name, centroids_files, args.num_cores)
            tsprint("MIDAS2::build_bowtie2_pangenomes_indexes::finish")

    except Exception as error:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error. Specify --debug flag to keep.")
            command(f"rm -f {args.bt2_indexes_dir}")
        raise error



@register_args
def main(args):
    tsprint(f"Build Bowtie2 genome database in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    build_bowtie2db(args)
