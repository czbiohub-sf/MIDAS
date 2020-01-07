import json
import os
from collections import defaultdict
import numpy as np
import Bio.SeqIO
from pysam import AlignmentFile  # pylint: disable=no-name-in-module

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, num_physical_cores, command, InputStream, OutputStream, select_from_tsv, multithreading_hashmap, download_reference, split
from iggtools.params import outputs


DEFAULT_ALN_COV = 0.75
DEFAULT_SPECIES_COVERAGE = 3.0
DEFAULT_ALN_MAPID = 94.0
DEFAULT_ALN_MAPQ = 20
DEFAULT_ALN_READQ = 20
DEFAULT_ALN_BASEQ = 30
DEFAULT_ALN_TRIM = 0

DECIMALS = ".6f"


def register_args(main_func):
    subparser = add_subcommand('midas_run_snps', main_func, help='single-nucleotide-polymorphism prediction')
    subparser.add_argument('outdir',
                           type=str,
                           help="""Path to directory to store results.  Name should correspond to unique sample identifier.""")
    subparser.add_argument('-1',
                           dest='r1',
                           required=True,
                           help="FASTA/FASTQ file containing 1st mate if using paired-end reads.  Otherwise FASTA/FASTQ containing unpaired reads.")
    subparser.add_argument('-2',
                           dest='r2',
                           help="FASTA/FASTQ file containing 2nd mate if using paired-end reads.")
    subparser.add_argument('--max_reads',
                           dest='max_reads',
                           type=int,
                           metavar="INT",
                           help=f"Number of reads to use from input file(s).  (All)")
    subparser.add_argument('--species_cov',
                           type=float,
                           dest='species_cov',
                           metavar='FLOAT',
                           default=DEFAULT_SPECIES_COVERAGE,
                           help=f"Include species with >X coverage ({DEFAULT_SPECIES_COVERAGE})")
    if False:
        # This is unused.
        subparser.add_argument('--species_topn',
                               type=int,
                               dest='species_topn',
                               metavar='INT',
                               help='Include top N most abundant species')

    # Alignment flags
    subparser.add_argument('--aln_speed',
                           type=str,
                           dest='aln_speed',
                           default='very-sensitive',
                           choices=['very-fast', 'fast', 'sensitive', 'very-sensitive'],
                           help='Alignment speed/sensitivity (very-sensitive).  If aln_mode is local (default) this automatically issues the corresponding very-sensitive-local, etc flag to bowtie2.')
    subparser.add_argument('--aln_mode',
                           type=str,
                           dest='aln_mode',
                           default='global',
                           choices=['local', 'global'],
                           help='Global/local read alignment (local, corresponds to the bowtie2 --local).  Global corresponds to the bowtie2 default --end-to-end.')
    subparser.add_argument('--aln_interleaved',
                           action='store_true',
                           default=False,
                           help='FASTA/FASTQ file in -1 are paired and contain forward AND reverse reads')

    #  Pileup flags (samtools, or postprocessing)
    subparser.add_argument('--aln_mapid',
                           dest='aln_mapid',
                           type=float,
                           metavar="FLOAT",
                           default=DEFAULT_ALN_MAPID,
                           help=f"Discard reads with alignment identity < MAPID.  Values between 0-100 accepted.  ({DEFAULT_ALN_MAPID})")
    subparser.add_argument('--aln_mapq',
                           dest='aln_mapq',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_ALN_MAPQ,
                           help=f"Discard reads with alignment identity < MAPQ. ({DEFAULT_ALN_MAPID})")
    subparser.add_argument('--aln_readq',
                           dest='aln_readq',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_ALN_READQ,
                           help=f"Discard reads with mean quality < READQ ({DEFAULT_ALN_READQ})")
    subparser.add_argument('--aln_cov',
                           dest='aln_cov',
                           default=DEFAULT_ALN_COV,
                           type=float,
                           metavar="FLOAT",
                           help=f"Discard reads with alignment coverage < ALN_COV ({DEFAULT_ALN_COV}).  Values between 0-1 accepted.")
    subparser.add_argument('--aln_baseq',
                           dest='aln_baseq',
                           default=DEFAULT_ALN_BASEQ,
                           type=int,
                           metavar="INT",
                           help=f"Discard bases with quality < ALN_BASEQ ({DEFAULT_ALN_BASEQ})")
    subparser.add_argument('--aln_trim',
                           dest='aln_trim',
                           default=DEFAULT_ALN_TRIM,
                           type=int,
                           metavar="INT",
                           help=f"Trim ALN_TRIM base-pairs from 3'right end of read ({DEFAULT_ALN_TRIM})")
    subparser.add_argument('--aln_discard',
                           dest='aln_baq',
                           default=False,
                           help='Discard discordant read-pairs (False)')
    subparser.add_argument('--aln_baq',
                           dest='aln_baq',
                           default=False,
                           help='Enable BAQ: per-base alignment quality (False)')
    subparser.add_argument('--aln_adjust_mq',
                           dest='aln_adjust_mq',
                           default=False,
                           help='Adjust MAPQ (False)')
    return main_func



def midas_run_snps(args):
    print("hello world")



@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_run_snps(args)
