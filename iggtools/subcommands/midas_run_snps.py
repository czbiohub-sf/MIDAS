import json
import os
from collections import defaultdict
import multiprocessing

import numpy as np
from pysam import AlignmentFile  # pylint: disable=no-name-in-module
import Bio.SeqIO

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, num_physical_cores, command, InputStream, OutputStream, multithreading_hashmap, download_reference
from iggtools.params import outputs
from iggtools.models.uhgg import UHGG
from iggtools.common.samples import parse_species_profile, select_species
from iggtools.common.bowtie2 import build_bowtie2_db, bowtie2_align, samtools_index


DEFAULT_ALN_COV = 0.75
DEFAULT_SPECIES_COVERAGE = 3.0
DEFAULT_ALN_MAPID = 94.0
DEFAULT_ALN_MAPQ = 20
DEFAULT_ALN_READQ = 20
DEFAULT_ALN_BASEQ = 30
DEFAULT_ALN_TRIM = 0

DECIMALS = ".3f"


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
                           help='Global/local read alignment (local, corresponds to the bowtie2 --local; global corresponds to the bowtie2 default --end-to-end).')
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
                           help=f"Discard reads with mapping quality < MAPQ. ({DEFAULT_ALN_MAPQ})")
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
    if False:
        # This is unused.
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
    subparser.add_argument('--sparse',
                           dest='sparse',
                           default=False,
                           help=f"Omit zero rows from output.")
    return main_func


def imported_genome_file(genome_id, species_id, component):
    return f"{outputs.cleaned_imports}/{species_id}/{genome_id}/{genome_id}.{component}"


def keep_read_worker(aln, args, aln_stats):
    aln_stats['aligned_reads'] += 1

    align_len = len(aln.query_alignment_sequence)
    query_len = aln.query_length

    # min pid
    if 100 * (align_len - dict(aln.tags)['NM']) / float(align_len) < args.aln_mapid:
        return False
    # min read quality
    if np.mean(aln.query_qualities) < args.aln_readq:
        return False
    # min map quality
    if aln.mapping_quality < args.aln_mapq:
        return False
    # min aln cov
    if align_len / float(query_len) < args.aln_cov:
        return False

    aln_stats['mapped_reads'] += 1
    return True


def species_pileup(species_id, args, tempdir, outputdir, contig_file, contigs_db_stats):
    # Read in contigs information for current species_id

    contigs = {}
    contigs_db_stats['species_counts'] += 1 # not being updated and passed as expected

    with InputStream(contig_file) as file:
        for rec in Bio.SeqIO.parse(file, 'fasta'):
            contigs[rec.id] = {
                "species_id": species_id,
                "contig_len": int(len(rec.seq)),
                "contig_seq": str(rec.seq),
            }
            contigs_db_stats['total_length'] += contigs[rec.id]["contig_len"]
            contigs_db_stats['total_seqs'] += 1

    # Summary statistics
    aln_stats = {
        "genome_length": 0,
        "total_depth": 0,
        "covered_bases": 0,
        "aligned_reads":0,
        "mapped_reads":0,
        }

    def keep_read(x):
        return keep_read_worker(x, args, aln_stats)

    header = ['ref_id', 'ref_pos', 'ref_allele', 'depth', 'count_a', 'count_c', 'count_g', 'count_t']
    path = f"{outputdir}/{species_id}.snps.lz4"

    with OutputStream(path) as file:

        file.write('\t'.join(header) + '\n')
        zero_rows_allowed = not args.sparse

        # Loop over alignment for current species's contigs
        with AlignmentFile(f"{tempdir}/repgenomes.bam") as bamfile:
            for contig_id in sorted(list(contigs.keys())): # why need to sort?
                contig = contigs[contig_id]
                counts = bamfile.count_coverage(
                    contig_id,
                    start=0,
                    end=contig["contig_len"],
                    quality_threshold=args.aln_baseq,
                    read_callback=keep_read)

                for ref_pos in range(0, contig["contig_len"]):
                    ref_allele = contig["contig_seq"][ref_pos]
                    depth = sum([counts[nt][ref_pos] for nt in range(4)])
                    count_a = counts[0][ref_pos]
                    count_c = counts[1][ref_pos]
                    count_g = counts[2][ref_pos]
                    count_t = counts[3][ref_pos]
                    values = [contig_id, ref_pos + 1, ref_allele, depth, count_a, count_c, count_g, count_t]

                    if depth > 0 or zero_rows_allowed:
                        file.write('\t'.join(str(val) for val in values) + '\n')

                    aln_stats['genome_length'] += 1
                    aln_stats['total_depth'] += depth
                    if depth > 0:
                        aln_stats['covered_bases'] += 1

    tsprint(json.dumps({species_id: aln_stats}, indent=4))
    return (species_id, {k: str(v) for k, v in aln_stats.items()})


def pysam_pileup(args, species_ids, tempdir, outputdir, contigs_files):
    "Counting alleles and run pileups per species in parallel"

    # Update alignment stats for species
    species_pileup_stats = defaultdict()
    contigs_db_stats = {'species_counts':0, 'total_seqs':0, 'total_length':0}

    mp = multiprocessing.Pool(num_physical_cores)
    argument_list = [(sp_id, args, tempdir, outputdir, contigs_files[sp_id], contigs_db_stats) for sp_id in species_ids]

    for species_id, aln_stats in mp.starmap(species_pileup, argument_list):
        sp_stats = {
            "genome_length": int(aln_stats['genome_length']),
            "covered_bases": int(aln_stats['covered_bases']),
            "total_depth": int(aln_stats['total_depth']),
            "aligned_reads": int(aln_stats['aligned_reads']),
            "mapped_reads": int(aln_stats['mapped_reads']),
            "fraction_covered": 0.0,
            "mean_coverage": 0.0,
        }

        if sp_stats["genome_length"] > 0:
            sp_stats["fraction_covered"] = format(sp_stats["covered_bases"] / sp_stats["genome_length"], DECIMALS)

        if sp_stats["covered_bases"] > 0:
            sp_stats["mean_coverage"] = format(sp_stats["total_depth"] / sp_stats["covered_bases"], DECIMALS)

        species_pileup_stats[species_id] = sp_stats

    tsprint(f"contigs_db_stats - total genomes: {contigs_db_stats['species_counts']}")
    tsprint(f"contigs_db_stats - total contigs: {contigs_db_stats['total_seqs']}")
    tsprint(f"contigs_db_stats - total base-pairs: {contigs_db_stats['total_length']}")

    return species_pileup_stats


def write_snps_summary(species_pileup_stats, outfile):
    """ Get summary of mapping statistics """
    header = ['species_id', 'genome_length', 'covered_bases', 'total_depth', 'aligned_reads', 'mapped_reads', 'fraction_covered', 'mean_coverage']
    with OutputStream(outfile) as file:
        file.write('\t'.join(header) + '\n')
        for species_id, species_aln in species_pileup_stats.items():
            values = list(species_aln.values())
            values.insert(0, species_id)
            file.write('\t'.join(map(str, values)) + '\n')


def midas_run_snps(args):

    tempdir = f"{args.outdir}/snps/temp_sc{args.species_cov}"
    if args.debug and os.path.exists(tempdir):
        tsprint(f"INFO:  Reusing existing temp data in {tempdir} according to --debug flag.")
    else:
        command(f"rm -rf {tempdir}")
        command(f"mkdir -p {tempdir}")

    outputdir = f"{args.outdir}/snps/output_sc{args.species_cov}"
    if not os.path.exists(outputdir):
        command(f"mkdir -p {outputdir}")

    try:
        # The full species profile must exist -- it is output by run_midas_species.
        # Restrict to species above requested coverage.
        full_species_profile = parse_species_profile(args.outdir)
        species_profile = select_species(full_species_profile, args.species_cov)

        local_toc = download_reference(outputs.genomes)
        db = UHGG(local_toc)
        representatives = db.representatives

        def download_contigs(species_id):
            return download_reference(imported_genome_file(representatives[species_id], species_id, "fna.lz4"), f"{tempdir}/{species_id}")

        # Download repgenome_id.fna for every species in the restricted species profile.
        contigs_files = multithreading_hashmap(download_contigs, species_profile.keys(), num_threads=20)

        # Use Bowtie2 to map reads to a representative genomes
        bt2_db_name = "repgenomes"
        build_bowtie2_db(tempdir, bt2_db_name, contigs_files)
        bowtie2_align(args, tempdir, bt2_db_name, sort_aln=True)

        # Use mpileup to identify SNPs
        samtools_index(args, tempdir, bt2_db_name)
        species_pileup_stats = pysam_pileup(args, list(species_profile.keys()), tempdir, outputdir, contigs_files)

        write_snps_summary(species_pileup_stats, f"{args.outdir}/snps/output_sc{args.species_cov}/summary.txt")

    except:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error. Specify --debug flag to keep.")
            command(f"rm -rf {tempdir}", check=False)
            command(f"rm -rf {outputdir}", check=False)


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_run_snps(args)
