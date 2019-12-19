import json

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint
from iggtools import params


DEFAULT_WORD_SIZE = 28
DEFAULT_ALN_COV = 0.75


def register_args(main_func):
    subparser = add_subcommand('midas_run_species', main_func, help='compute species abundance profile for given sample(s)')
    subparser.add_argument('-1',
                           dest='r1',
                           required=True,
                           help="FASTA/FASTQ file containing 1st mate if using paired-end reads.  Otherwise FASTA/FASTQ containing unpaired reads.")
    subparser.add_argument('-2',
                           dest='r2',
                           help="FASTA/FASTQ file containing 2nd mate if using paired-end reads.")
    subparser.add_argument('--word_size',
                           dest='word_size',
                           default=DEFAULT_WORD_SIZE,
                           type=int,
                           metavar="INT",
                           help=f"Word size for BLAST search ({DEFAULT_WORD_SIZE}).  Use word sizes > 16 for greatest efficiency.")
    subparser.add_argument('--mapid',
                           dest='mapid',
                           type=float,
                           metavar="FLOAT",
                           help=f"Discard reads with alignment identity < MAPID.  Values between 0-100 accepted.  By default gene-specific species-level cutoffs are used, as specifeid in {params.inputs.marker_genes_hmm_cutoffs}")
    subparser.add_argument('--aln_cov',
                           dest='aln_cov',
                           default=DEFAULT_ALN_COV,
                           type=float,
                           metavar="FLOAT",
                           help=f"Discard reads with alignment coverage < ALN_COV ({DEFAULT_ALN_COV}).  Values between 0-1 accepted.")
    subparser.add_argument('--read_length',
                           dest='read_length',
                           type=int,
                           metavar="INT",
                           help=f"Trim reads to READ_LENGTH and discard reads with length < READ_LENGTH.  By default, reads are not trimmed or filtered")
    return main_func


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
