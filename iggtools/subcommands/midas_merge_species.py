import json
import os
from collections import defaultdict
import numpy as np

from iggtools.models.pool import Pool, select_species

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, command, InputStream, OutputStream, select_from_tsv
from iggtools.params.schemas import species_profile_schema, DEFAULT_SAMPLE_DEPTH, DECIMALS


def register_args(main_func):
    subparser = add_subcommand('midas_merge_species', main_func, help='merge MIDAS species abundance results across metagenomic samples')
    subparser.add_argument('outdir',
                           type=str,
                           help="""Path to directory to store results.  Name should correspond to unique sample identifier.""")
    subparser.add_argument('--sample_list',
                           dest='sample_list',
                           type=str,
                           required=True,
                           help=f"TSV file mapping sample name to midas_run_species.py output directories")
    subparser.add_argument('--sample_depth',
                           dest='sample_depth',  # min_cov
                           type=float,
                           metavar="FLOAT",
                           default=DEFAULT_SAMPLE_DEPTH,
                           help=f"Minimum per-sample marker-gene-depth for estimating species prevalence ({DEFAULT_SAMPLE_DEPTH})")
    if False:
        # This is not currently in use.
        subparser.add_argument('-input',
                               dest='input',
                               type=str,
                               required=True,
                               help=f"Input to sample directories output by run_midas.py; see '-t' for details")
        subparser.add_argument('-input_type',
                               dest='intype',
                               type=str,
                               required=True,
                               metavar="INPUT_TYPE",
                               choices=['list', 'file', 'dir'],
                               help=f"Specify one of the following: (list, dir, file) \n list: -i is a comma-separated list (eg: /samples/sample_1,/samples/sample_2) \n dir: -i is a directory containing all samples (eg: /samples) \n file: -i is a file of paths to samples (eg: /sample_paths.txt)")
        subparser.add_argument('--max_samples',
                               dest='max_samples',
                               type=int,
                               metavar="INT",
                               help=f"Maximum number of samples to process.Useful for testing (use all)")
    return main_func


def read_samples(sample_list):
    with InputStream(sample_list) as stream:
        samples = dict(select_from_tsv(stream, selected_columns=["sample_name", "midas_output_path"]))

    for sample_name in samples.keys():
        midas_output_path = samples[sample_name]
        assert os.path.exists(midas_output_path), f"MIDAS output directory {midas_output_path} for sample {sample_name} not exist."
        species_profile = f"{midas_output_path}/species/species_profile.txt"
        assert os.path.exists(species_profile), f"Missing MIDAS species profile: {species_profile}"
        samples[sample_name] = species_profile
    return samples


def read_species_profile(species_profile):
    with InputStream(species_profile) as instream:
        for record in select_from_tsv(instream, selected_columns=species_profile_schema, result_structure=dict):
            yield record


def transpose(samples, columns):
    """
    Collate given columns across samples and transpose the matrix by species_id
    samples: dict of sample_name: path_to_table_by_sample
    """
    transposed = defaultdict(dict)
    sample_counts = len(samples)

    for sample_index, species_profile in enumerate(samples.values()):
        for record in read_species_profile(species_profile):
            species_id = record["species_id"]
            for col in columns:
                acc = transposed[col].get(species_id)
                if not acc:
                    acc = [species_id] + [0.0] * sample_counts
                    transposed[col][species_id] = acc
                acc[1 + sample_index] = record[col]
    return transposed


def cal_prevalence(rowvector, threshold):
    return sum(1 if val >= threshold else 0 for val in rowvector)


def format_data(x):
    return format(x, DECIMALS) if isinstance(x, float) else str(x)


def compute_stats(tabundance, tcoverage, args):

    assert tabundance.keys() == tcoverage.keys(), f"Transposed abundance matrix and coverage matrix have different species ids"
    species_ids = tabundance.keys()

    stats = defaultdict()
    for species_id in species_ids:
        species_abun = tabundance[species_id][1:]
        species_cov = tcoverage[species_id][1:]
        values = {
            "species_id": species_id,
            "median_abundance": np.median(species_abun),
            "mean_abundance": np.mean(species_abun),
            "median_coverage": np.median(species_cov),
            "mean_coverage": np.mean(species_cov),
            "prevalence": cal_prevalence(species_cov, args.sample_depth)
        }
        stats[species_id] = values
    return stats


def write_results(samples, transposed, stats, outdir, sort_by="median_coverage"):
    "Write the transposed tables into separate files"

    header = ["species_id"] + list(samples.keys())
    for field in list(species_profile_schema.keys())[1:]:
        path = f"{outdir}/{field}.tsv"
        data = transposed[field]
        with OutputStream(path) as outfile:
            outfile.write("\t".join(header) + "\n")
            for values in data.values():
                outfile.write("\t".join(map(format_data, values)) + "\n")

    # Order species in stats by descending relative abundance
    outfile = f"{outdir}/species_prevalence.tsv"
    with OutputStream(outfile) as ostream:
        header = ["species_id", "median_abundance", "mean_abundance", "median_coverage", "mean_coverage", "prevalence"]
        ostream.write("\t".join(header) + "\n")

        sorted_species = sorted(((row[sort_by], row["species_id"]) for row in stats.values()), reverse=True)
        for spec_tuple in sorted_species:
            species_id = spec_tuple[1]
            values = map(format_data, list(stats[species_id].values()))
            ostream.write("\t".join(values) + "\n")


def midas_merge_species(args):

    outdir = f"{args.outdir}/merged/species"
    command(f"rm -rf {outdir}")
    command(f"mkdir -p {outdir}")

    pool_of_samples = Pool(args.sample_list)

    exit(0)
    # Slice the across-samples species profile matrix by species_id
    transposed = transpose(samples, list(species_profile_schema.keys())[1:])

    # Calculate summary statistics for coverage and relative abundance
    stats = compute_stats(transposed["relative_abundance"], transposed["coverage"], args)

    write_results(samples, transposed, stats, outdir, "median_coverage")



@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_merge_species(args)
