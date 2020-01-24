import json
import random
import os
from collections import defaultdict
from itertools import chain
import numpy as np

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, num_physical_cores, command, InputStream, OutputStream, select_from_tsv, multithreading_map, download_reference, TimedSection
from iggtools.models.uhgg import UHGG
from iggtools import params


DEFAULT_SAMPLE_DEPTH = 1.0


species_abundance_schema = {
    "species_id": str,
    "count_reads": int,
    "coverage": float,
    "relative_abundance": float
}


def register_args(main_func):
    subparser = add_subcommand('midas_merge_species', main_func, help='merge MIDAS species abundance results across metagenomic samples')
    subparser.add_argument('outdir',
                           type=str,
                           help="""Path to directory to store results.  Name should correspond to unique sample identifier.""")
    if False:
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
        subparser.add_argument('--max_samples',
                               dest='max_samples',
                               type=int,
                               metavar="INT",
                               help=f"Maximum number of samples to process.Useful for testing (use all)")
    return main_func



def decode_indirs_arg(args):

    assert args.intype in ["dir", "file", "list"], f"Specified input type {args.intype} is not valid."
    indirs = []

    if args.intype == "dir":
        assert os.path.isdir(args.input), f"Specied input directory {args.input} does not exist."
        for indir in os.listdir(args.input):
            assert os.path.isdir(indir), f"Specified input file {indir} does not exist."
            indirs.append(f"{args.input}/{indir}")

    if args.intype == "file":
        assert os.path.isfile(args.input), f"Specified input file {args.input} does not exist."
        with InputStream(args.input) as fin:
            for line in fin:
                indir = line.rstrip().rstrip('/')
                assert os.path.isdir(indir), f"Specified input file {indir} does not exist."
                indirs.append(indir)

    if args.intype == "list":
        for indir in args.input.split(','):
            assert os.path.isdir(indir), f"Specified input file {indir} does not exist."
            indirs.append(indir)

    return indirs


def identify_samples(sample_list):

    with InputStream(sample_list) as stream:
        samples = dict(select_from_tsv(stream, selected_columns={"sample_name":str, "midas_output_path":str}))

    for sample_name in samples.keys():
        midas_output_path = samples[sample_name]
        assert os.path.exists(midas_output_path), f"MIDAS output directory {midas_output_path} for sample {sample_name} not exist."

        species_profile = f"{midas_output_path}/species/species_profile.txt"
        assert os.path.exists(species_profile), f"Missing MIDAS species profile: {species_profile}"

        samples["species_profile"] = species_profile

    return samples


def read_abundance(species_profile_dir):
    "Return map of species_id to coverage for the species present in the sample."
    with InputStream(f"{species_profile_dir}") as stream:
        return dict(select_from_tsv(stream, selected_columns=species_abundance_schema))


def store_data(args, samples, species_ids):
    # Initialize data of dict(species_id) = dict()
    data = []
    for species_id in species_ids:
        data[species_id] = {}
        for field in ['relative_abundance', 'coverage', 'count_reads']:
            data[species_id][field] = []

    # I think this one liner should be able to
    data = defaultdict(lambda: defaultdict(list))

    # Read species_profile results for ALL specified samples for ALL the species
    for sample_name, species_profile in samples.items():
        # samples.parse_species_profile()
        abundance = read_abundance(species_profile)
        for species_id, values in abundance.items():
            for field in ['relative_abundance', 'coverage', 'count_reads']:
                if field in values:
                    data[species_id][field].append(values[field])
    # data[species_id] has three separate tablles
    # my question is: how do we keep track of the samples information
    return data


def prevalence(x, y):
    return sum(1 if rc >= y else 0 for rc in x)


def compute_stats(args, data):
    species_ids = data.keys()
    stats = dict((species_id, {}) for species_id in species_ids)

    for species_id in species_ids:
        # abundance
        x = data[species_id]['relative_abundance']
        stats[species_id]["median_abundance"] = np.median(x)
        stats[species_id]["mean_abundance"] = np.mean(x)
        # coverage
        x = data[species_id]['coverage']
        stats[species_id]["median_coverage"] = np.median(x)
        stats[species_id]["mean_coverage"] = np.mean(x)
        # prevalence
        stats[species_id]["prevalence"] = prevalence(x, args.sample_depth) #min_cov
    return stats



def write_abundance(args, samples, data):
    for field in ['relative_abundance', 'coverage', 'count_reads']:
        with OutputStream(f"{args.outdir}/{field}.txt") as outfile:
            outfile.write("\t".join(["species_id"] + list(samples.keys())) + "\n")
            for species_id in data:
                outfile.write(species_id)
                for x in data[species_id][field]:
                    outfile.write("\t%s" % str(x))
                outfile.write("\n")
    ## STOP here: the sample information is coded in the order of the data[species_id]['counts'] as the same order of the sample_names.
    ## todo: what if species_a doesn't exist in sample_b ?? !!

def midas_merge_species(args):
    # List samples and species

    #sample_dirs = decode_indirs_arg(args)
    samples = identify_samples(args.sample_list)

    db = UHGG()
    species_info = db.species

    # Read in data and computate stats






@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_merge_species(args)
