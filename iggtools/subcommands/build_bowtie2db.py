import json

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, num_physical_cores, command
from iggtools.common.bowtie2 import build_bowtie2_db
from iggtools.models.midasdb import MIDAS_DB
from iggtools.models.species import parse_species, filter_species


def register_args(main_func):
    subparser = add_subcommand('build_bowtie2db', main_func, help='build repgenome and pangenome bowtie2 indexes given list of species')

    subparser.add_argument('--midas_db',
                           dest='midas_db',
                           type=str,
                           metavar="CHAR",
                           required=True,
                           help=f"Local MIDAS DB which mirrors the s3 IGG db")
    subparser.add_argument('--bt2_indexes_dir',
                           dest='bt2_indexes_dir',
                           type=str,
                           metavar="CHAR",
                           required=True,
                           help=f"Path to bowtie2 indexes directory")
    subparser.add_argument('--bt2_indexes_name',
                           dest='bt2_indexes_name',
                           type=str,
                           metavar="CHAR",
                           default="repgenomes",
                           choices=['repgenomes', 'pangenomes'],
                           help=f"Bowtie2 db name to build")

    subparser.add_argument('--species_list',
                           dest='species_list',
                           type=str,
                           metavar="CHAR",
                           help=f"Comma separated list of species ids")

    subparser.add_argument('--species_profile',
                           dest='species_profile',
                           type=str,
                           metavar="CHAR",
                           help=f"Path to species coverage TSV file with species_id")
    subparser.add_argument('--select_by',
                           dest='select_by',
                           type=str,
                           metavar="CHAR",
                           default="sample_counts",
                           choices=['sample_counts', 'median_marker_coverage'],
                           help=f"Column from species_profile based on which to select species.")
    subparser.add_argument('--select_threshold',
                           dest='select_threshold',
                           type=float,
                           metavar="FLOAT",
                           help=f"Minimum threshold values of for selected columns.")

    subparser.add_argument('--num_cores',
                           dest='num_cores',
                           type=int,
                           metavar="INT",
                           default=num_physical_cores,
                           help=f"Number of physical cores to use ({num_physical_cores})")
    return main_func


def build_bowtie2db(args):

    try:
        if args.species_list:
            species_ids_of_interest = parse_species(args)
        elif args.species_profile and args.select_by and args.select_threshold:
            species_ids_of_interest = filter_species(args.species_profile, args.select_by, args.select_threshold)
        else:
            raise Exception(f"Need to provide either species_list or species_profile as input arguments")
        tsprint(f"CZ::build_bowtie2db::build bt2 indexees for the listed species: {species_ids_of_interest}")

        # Fetch UHGG related files
        midas_db = MIDAS_DB(args.midas_db, args.num_cores)

        if args.bt2_indexes_name == "repgenomes":
            tsprint(f"CZ::build_bowtie2_repgenomes_indexes::start")
            contigs_files = midas_db.fetch_files("representative_genome", species_ids_of_interest)
            tsprint(contigs_files)
            build_bowtie2_db(args.bt2_indexes_dir, args.bt2_indexes_name, contigs_files, args.num_cores)
            tsprint(f"CZ::build_bowtie2_repgenomes_indexes::finish")

        if args.bt2_indexes_name == "pangenomes":
            tsprint(f"CZ::build_bowtie2_pangenomes_indexes::start")
            centroids_files = midas_db.fetch_files("pangenome_centroids", species_ids_of_interest)
            tsprint(centroids_files)
            build_bowtie2_db(args.bt2_indexes_dir, args.bt2_indexes_name, centroids_files, args.num_cores)
            tsprint(f"CZ::build_bowtie2_pangenomes_indexes::finish")

    except Exception as error:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error. Specify --debug flag to keep.")
            command(f"rm -f {args.bt2_indexes_dir}")
        raise error


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    build_bowtie2db(args)
