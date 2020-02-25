import json
import os
from collections import defaultdict
import numpy as np

from iggtools.models.pool import Pool, select_species

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, command, InputStream, OutputStream, select_from_tsv
from iggtools.params.schemas import species_profile_schema, species_prevalence_schema, DECIMALS, fetch_default_genome_depth


DEFAULT_GENOME_DEPTH = fetch_default_genome_depth("species")

def register_args(main_func):
    subparser = add_subcommand('midas_merge_species', main_func, help='merge MIDAS species abundance results across metagenomic samples')
    subparser.add_argument('outdir',
                           type=str,
                           help="""Path to directory to store results.  Name should correspond to unique sample identifier.""")
    subparser.add_argument('--samples_list',
                           dest='samples_list',
                           type=str,
                           required=True,
                           help=f"TSV file mapping sample name to midas_run_species.py output directories")
    subparser.add_argument('--genome_depth',
                           dest='genome_depth',
                           type=float,
                           metavar="FLOAT",
                           default=DEFAULT_GENOME_DEPTH,
                           help=f"Minimum per-sample marker-gene-depth for estimating species prevalence ({DEFAULT_GENOME_DEPTH})")
    return main_func


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


def compute_prevalence(rowvector, threshold):
    return sum(1 if val >= threshold else 0 for val in rowvector)


def format_data(x):
    return format(x, DECIMALS) if isinstance(x, float) else str(x)


def compute_stats(tabundance, tcoverage, args):

    assert tabundance.keys() == tcoverage.keys(), f"midas_merge_species::compute_stats() merged abun and cov matrices have different species_id orders"
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


def write_results(pool_of_samples, transposed, stats, sort_by="median_coverage"):
    """ Write the transposed tables into separate files """

    sample_names = pool_of_samples.fetch_sample_names()
    output_files = pool_of_samples.fetch_merged_species_output_files()

    for col, path in output_files.items():
        data = transposed[col]
        with OutputStream(path) as outfile:
            outfile.write("\t".join(["species_id"] + sample_names) + "\n")
            for values in data.values():
                outfile.write("\t".join(map(format_data, values)) + "\n")

    # Sort species in stats by descending relative abundance
    with OutputStream(pool_of_samples.fetch_species_prevalence_path()) as ostream:
        ostream.write("\t".join(list(species_prevalence_schema.keys())) + "\n")

        c_sort_by = list(species_prevalence_schema.keys()).index(sort_by)
        sorted_species = sorted(((row[c_sort_by], species_id) for species_id, row in stats.items()), reverse = True)

        for species_tuple in sorted_species:
            species_id = species_tuple[1]
            values = map(format_data, stats[species_id])
            ostream.write("\t".join(values) + "\n")


def midas_merge_species(args):

    paramstr = f"gd{args.genome_depth}"
    pool_of_samples = Pool(args.samples_list, args.outdir, paramstr, "species")

    # Slice the across-samples species profile matrix by species_id
    transposed = transpose(pool_of_samples, list(species_profile_schema.keys())[1:])

    # Calculate summary statistics for coverage and relative abundance
    stats = compute_stats(transposed["relative_abundance"], transposed["coverage"], args)
    write_results(pool_of_samples, transposed, stats, "median_coverage")



@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_merge_species(args)
