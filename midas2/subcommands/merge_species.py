#!/usr/bin/env python3
import json
from collections import defaultdict
import numpy as np

from midas2.models.samplepool import SamplePool
from midas2.common.argparser import add_subcommand
from midas2.common.utils import tsprint, OutputStream
from midas2.params.schemas import species_merge_schema, species_prevalence_schema, format_data


DEFAULT_MARKER_DEPTH = 0.0


def register_args(main_func):
    subparser = add_subcommand('merge_species', main_func, help='merge MIDAS species abundance results across metagenomic samples')
    subparser.add_argument('midas_outdir',
                           type=str,
                           help="""Path to directory to store results.  Name should correspond to unique sample identifier.""")
    subparser.add_argument('--samples_list',
                           dest='samples_list',
                           type=str,
                           required=True,
                           help=f"TSV file mapping sample name to run_species.py output directories")
    subparser.add_argument('--marker_depth',
                           dest='marker_depth',
                           type=float,
                           metavar="FLOAT",
                           default=DEFAULT_MARKER_DEPTH,
                           help=f"Minimum per-sample marker-gene-depth for estimating species prevalence ({DEFAULT_MARKER_DEPTH})")
    return main_func


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
              compute_prevalence(species_cov, args.marker_depth)]

        stats[species_id] = values
    return stats


def write_stats(stats, species_prevalence_filepath, sort_by="median_coverage"):
    """ Sort species in stats by descending relative abundance """
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
    col_names = list(species_merge_schema.keys())[1:]

    for col in col_names:
        outpath = pool_of_samples.get_target_layout(f"species_{col}")
        with OutputStream(outpath) as outfile:
            outfile.write("\t".join(["species_id"] + sample_names) + "\n")
            for values in transposed[col].values():
                outfile.write("\t".join(map(format_data, values)) + "\n")


def merge_species(args):

    try:
        global global_args
        global_args = args

        # Load species_summary into sample.profile
        pool_of_samples = SamplePool(args.samples_list, args.midas_outdir, "species")
        pool_of_samples.init_samples("species")
        pool_of_samples.create_dirs(["outdir"], args.debug)

        # Slice the across-samples species profile matrix by species_id
        tsprint(f"MIDAS2::write_species_results::start")
        cols = list(species_merge_schema.keys())[1:]
        transposed = transpose(pool_of_samples, cols)
        write_species_results(pool_of_samples, transposed)
        tsprint(f"MIDAS2::write_species_results::finish")

        # Calculate summary statistics for coverage and relative abundance
        tsprint(f"MIDAS2::write_stats::start")
        stats = compute_stats(transposed["marker_relative_abundance"], transposed["median_marker_coverage"]) #<-- marker_coverage
        write_stats(stats, pool_of_samples.get_target_layout("species_prevalence"), "median_coverage")
        tsprint(f"MIDAS2::write_stats::finish")

    except Exception as error:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error. Specify --debug flag to keep.")
            pool_of_samples.remove_dirs(["outdir"])
        raise error


@register_args
def main(args):
    tsprint(f"Merge species coverage profile in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    merge_species(args)
