import json
import os
from collections import defaultdict
import numpy as np

from iggtools.models.pool import Pool
from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, command, InputStream, OutputStream, select_from_tsv, command, download_reference
from iggtools.params.schemas import species_profile_schema, species_prevalence_schema, DECIMALS, fetch_default_genome_depth, format_data
from iggtools.params import outputs
from iggtools.common.bowtie2 import build_bowtie2_db, bowtie2_align, samtools_index, bowtie2_index_exists
from iggtools.models.uhgg import UHGG


DEFAULT_GENOME_DEPTH = fetch_default_genome_depth("species")


def compute_prevalence(rowvector, threshold):
    return sum(1 if val >= threshold else 0 for val in rowvector)


def transpose(pool_of_samples, columns):
    """ Collect given columns across samples and transpose the matrix by species_id """
    transposed = defaultdict(dict)

    total_samples_count = len(pool_of_samples.samples)

    for sample_index, sample in enumerate(pool_of_samples.samples):
        for species_id, species_record in sample.profile.items():
            for col in columns:
                acc = transposed[col].get(species_id)
                if not acc:
                    acc = [species_id] + [0.0] * total_samples_count
                    transposed[col][species_id] = acc
                acc[1 + sample_index] = species_record[col]
    return transposed


def compute_stats(tabundance, tcoverage):
    global global_args
    args = global_args

    assert tabundance.keys() == tcoverage.keys(), f"compute_and_write_stats::merged abun and cov matrices have different species_id orders"
    species_ids = tabundance.keys()
    stats = defaultdict()
    for species_id in species_ids:
        species_abun = tabundance[species_id][1:]
        species_cov = tcoverage[species_id][1:]

        values = [species_id, np.median(species_abun), np.mean(species_abun), \
              np.median(species_cov), np.mean(species_cov), \
              compute_prevalence(species_cov, args.genome_depth)]

        stats[species_id] = values
    return stats


def write_stats(stats, species_prevalence_filepath, sort_by="median_coverage"):
    # Sort species in stats by descending relative abundance
    with OutputStream(species_prevalence_filepath) as ostream:

        colnames = list(species_prevalence_schema.keys())
        ostream.write("\t".join(colnames) + "\n")

        c_sort_by = colnames.index(sort_by)
        sorted_species = sorted(((row[c_sort_by], species_id) for species_id, row in stats.items()), reverse=True)

        for species_tuple in sorted_species:
            species_id = species_tuple[1]
            ostream.write("\t".join(map(format_data, stats[species_id])) + "\n")


def write_species_results(pool_of_samples, transposed):
    """ Write the transposed tables into separate files """

    sample_names = pool_of_samples.fetch_samples_names()

    col_names = list(species_profile_schema.keys())[1:]

    for col in col_names:
        outpath = pool_of_samples.get_target_layout(f"species_{col}")

        with OutputStream(outpath) as outfile:
            outfile.write("\t".join(["species_id"] + sample_names) + "\n")
            for values in transposed[col].values():
                outfile.write("\t".join(map(format_data, values)) + "\n")


def midas_merge_species(args):

    try:
        global global_args
        global_args = args

        pool_of_samples = Pool(args.samples_list, args.midas_outdir, "species")
        # load species_summary into sample.profile
        pool_of_samples.init_samples("species")
        # create output and temp directory
        pool_of_samples.create_dirs(["outdir", "tempdir"], args.debug)

        # Slice the across-samples species profile matrix by species_id
        cols = list(species_profile_schema.keys())[1:]
        transposed = transpose(pool_of_samples, cols)
        write_species_results(pool_of_samples, transposed)

        # Calculate summary statistics for coverage and relative abundance
        stats = compute_stats(transposed["rel_abundance"], transposed["coverage"])
        write_stats(stats, pool_of_samples.get_target_layout("species_prevalence"), "median_coverage")

        # TO move to another subcommand
        if args.build_bowtie2_db:
            # The input for this section is species_prevalance.tsv
            species_ids_of_interest = []
            for species_id, record in stats.items():
                if record[-1] > 0:
                    species_ids_of_interest.append(species_id)

            # Create the dbs/species
            pool_of_samples.create_species_subdir(species_ids_of_interest, "dbsdir", args.debug)
            pool_of_samples.create_species_subdir(species_ids_of_interest, "dbs_tempdir", args.debug)

            # where to build the shared dbs for all the samples to merge
            bt2_db_dir = pool_of_samples.get_target_layout("dbsdir")
            bt2_db_temp_dir = pool_of_samples.get_target_layout("dbs_tempdir")

            rep_bt2_db_name = "repgenomes"
            pan_bt2_db_name = "pangenomes"

            #
            local_toc = download_reference(outputs.genomes, bt2_db_temp_dir)
            db = UHGG(local_toc)

            # Fetch the files per genomes
            contigs_files = db.fetch_files(species_ids_of_interest, bt2_db_temp_dir, filetype="contigs")
            centroids_files = db.fetch_files(species_ids_of_interest, bt2_db_temp_dir, filetype="centroids")

            if False:
                build_bowtie2_db(bt2_db_dir, rep_bt2_db_name, contigs_files)
                build_bowtie2_db(bt2_db_dir, pan_bt2_db_name, centroids_files)

            # What do we need and what is actually we are doing?
    except:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error. Specify --debug flag to keep.")
            pool_of_samples.remove_dirs(["outdir", "tempdir", "dbsdir"])
        raise



def register_args(main_func):
    subparser = add_subcommand('midas_merge_species', main_func, help='merge MIDAS species abundance results across metagenomic samples')
    subparser.add_argument('midas_outdir',
                           type=str,
                           help="""Path to directory to store results.  Name should correspond to unique sample identifier.""")
    subparser.add_argument('--samples_list',
                           dest='samples_list',
                           type=str,
                           required=True,
                           help=f"TSV file mapping sample name to midas_run_species.py output directories")
    subparser.add_argument('--species_list',
                           dest='species_list',
                           type=str,
                           metavar="CHAR",
                           help=f"Comma separated list of species ids")
    subparser.add_argument('--genome_depth',
                           dest='genome_depth',
                           type=float,
                           metavar="FLOAT",
                           default=DEFAULT_GENOME_DEPTH,
                           help=f"Minimum per-sample marker-gene-depth for estimating species prevalence ({DEFAULT_GENOME_DEPTH})")
    subparser.add_argument('--build_bowtie2_db',
                           action='store_true',
                           default=False,
                           help=f"Omit zero rows from output.")
    return main_func


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_merge_species(args)
