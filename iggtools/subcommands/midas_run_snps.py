import json
import os
from collections import defaultdict
import multiprocessing

import numpy as np
from pysam import AlignmentFile  # pylint: disable=no-name-in-module
import Bio.SeqIO

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, num_physical_cores, command, InputStream, OutputStream, select_from_tsv, multithreading_hashmap, download_reference, split
from iggtools.params import outputs
from iggtools.models.uhgg import UHGG


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
    subparser.add_argument('--sparse',
                           dest='sparse',
                           default=False,
                           help=f"Omit zero rows from output.")
    return main_func


def parse_species_profile(outdir):
    "Return map of species_id to coverage for the species present in the sample."
    with InputStream(f"{outdir}/species/species_profile.txt") as stream:
        return dict(select_from_tsv(stream, {"species_id": str, "coverage": float}))


def select_species(species_profile, coverage_threshold):
    return {species_id: species_coverage for species_id, species_coverage in species_profile.items() if species_coverage >= coverage_threshold}


def imported_genome_file(genome_id, species_id, component):
    return f"{outputs.cleaned_imports}/{species_id}/{genome_id}/{genome_id}.{component}"


def initialize_contigs(contigs_files):
    contigs = {}
    db_stats = {'total_length':0, 'total_seqs':0, 'species':0}
    for species_id, contig_filename in contigs_files.items():
        db_stats['species'] += 1
        with InputStream(contig_filename) as file:
            for rec in Bio.SeqIO.parse(file, 'fasta'):
                contigs[rec.id] = {
                    "contig_seq": str(rec.seq),
                    "contig_len": len(rec.seq),
                    "species_id": species_id,
                }
                db_stats['total_length'] += contigs[rec.id]["contig_len"]
                db_stats['total_seqs'] += 1
    return contigs, db_stats


def build_repgenome_db(tempdir, contigs):
    "Build Bowtie2 database of representative genomes for the species present in the sample"
    if all(os.path.exists(f"{tempdir}/repgenomes.{ext}") for ext in ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]):
        tsprint("Skipping bowtie2-build as database files appear to exist.")
        return
    command(f"rm -f {tempdir}/repgenomes.fa")
    command(f"touch {tempdir}/repgenomes.fa")
    for files in split(contigs.values(), 20):  # keep "cat" commands short
        command("cat " + " ".join(files) + f" >> {tempdir}/repgenomes.fa")
    command(f"bowtie2-build --threads {num_physical_cores} {tempdir}/repgenomes.fa {tempdir}/repgenomes > {tempdir}/bowtie2-build.log")


def repgenome_align(args, tempdir):
    "Use Bowtie2 to map reads to specified representative genomes"
    if args.debug and os.path.exists(f"{tempdir}/repgenomes.bam"):
        tsprint(f"Skipping alignment in debug mode as temporary data exists: {tempdir}/repgenomes.bam")
        return

    max_reads = f"-u {args.max_reads}" if args.max_reads else ""
    aln_mode = "local" if args.aln_mode == "local" else "end-to-end"
    aln_speed = args.aln_speed if aln_mode == "end_to_end" else args.aln_speed + "-local"
    r2 = ""
    if args.r2:
        r1 = f"-1 {args.r1}"
        r2 = f"-2 {args.r2}"
    elif args.aln_interleaved:
        r1 = f"--interleaved {args.r1}"
    else:
        r1 = f"-U {args.r1}"

    try:
        command(f"set -o pipefail; bowtie2 --no-unal -x {tempdir}/repgenomes {max_reads} --{aln_mode} --{aln_speed} --threads {num_physical_cores} -q {r1} {r2} | \
                samtools view --threads {num_physical_cores} -b - | \
                samtools sort --threads {num_physical_cores} -o {tempdir}/repgenomes.bam")
    except:
        tsprint("Repgnome align run into error")
        command(f"rm -f {tempdir}/repgenomes.bam")
        raise


def samtools_index(tempdir, args):
    if args.debug and os.path.exists(f"{tempdir}/repgenomes.bam.bai"):
        tsprint(f"Skipping alignment in debug mode as temporary data exists: {tempdir}/repgenomes.bam")
        return
    try:
        command(f"samtools index -@ {num_physical_cores} {tempdir}/repgenomes.bam")
    except:
        command(f"rm -f {tempdir}/repgenomes.bam.bai")
        raise


def keep_read_worker(aln, min_pid, min_readq, min_mapq, min_aln_cov, aln_stats):
    aln_stats['aligned_reads'] += 1

    align_len = len(aln.query_alignment_sequence)
    query_len = aln.query_length
    # min pid
    if 100 * (align_len - dict(aln.tags)['NM']) / float(align_len) < min_pid:
        return False
    # min read quality
    if np.mean(aln.query_qualities) < min_readq:
        return False
    # min map quality
    if aln.mapping_quality < min_mapq:
        return False
    # min aln cov
    if align_len / float(query_len) < min_aln_cov:
        return False
    return True


def species_pileup(species_id):

    global global_args
    args = global_args

    global global_contigs
    contigs = global_contigs

    # Summary statistics
    aln_stats = {
        "genome_length": 0,
        "total_depth": 0,
        "covered_bases": 0,
        "aligned_reads":0,
        "mapped_reads":0,
        }

    def keep_read(x):
        return keep_read_worker(x, args.min_pid, args.min_readq, args.min_mapq, args.min_aln_cov, aln_stats)

    path = f"{args.outdir}/snps/output/{species_id}.snps.lz4"
    header = ['ref_id', 'ref_pos', 'ref_allele', 'depth', 'count_a', 'count_c', 'count_g', 'count_t']
    with OutputStream(path) as file:
        file.write('\t'.join(header) + '\n')

        zero_rows_allowed = not args['sparse']

        # loog over alignments for current species_id
        tempdir = f"{args.outdir}/snps/temp_sc{args.species_cov}" # idealy should pass on as parameter
        with AlignmentFile(f"{tempdir}/repgenomes.bam") as bamfile:
            for contig_id in sorted(list(contigs.keys())): # why do we need to sorted ?

                if contigs[contig_id]["species_id"] != species_id:
                    continue

                # Unpack frequently used variables
                contig_length = int(contigs[contig_id]["contig_length"])
                contig_seq = contigs[contig_id]["contig_seq"]

                counts = bamfile.count_coverage(
                    contig_id,
                    start=0,
                    end=contig_length,
                    quality_threshold=args.aln_baseq,
                    read_callback=keep_read)

                for i in range(0, contig_length):
                    ref_pos = i + 1
                    ref_allele = contig_seq[i]
                    depth = sum([counts[_][i] for _ in range(4)])
                    count_a = counts[0][i]
                    count_c = counts[1][i]
                    count_g = counts[2][i]
                    count_t = counts[3][i]
                    values = [contig_id, ref_pos, ref_allele, depth, count_a, count_c, count_g, count_t]
                    if depth > 0 or zero_rows_allowed:
                        file.write('\t'.join([str(val) for val in values])+'\n')
                    aln_stats['genome_length'] += 1
                    aln_stats['total_depth'] += depth
                    if depth > 0:
                        aln_stats['covered_bases'] += 1

    return (species_id, {k: str(v) for k, v in aln_stats.items()})


def pysam_pileup(args, species_ids, contigs):
    # Counting alleles

    # We cannot pass args to a subprocess unfortunately because args['log'] is an object;
    # so we can make it a global, although that is certainly living dangerously.
    # Can we still not be able to pass args as an object??

    global global_args
    global_args = args

    global global_contigs
    global_contigs = contigs

    # Run pileups per species in parallel
    # We might not need this for contigs.  It was an attempt to eliminate the nonserializable subprocess argument.  Which is args.

    # Update alignment stats for species
    species_alnstats = defaultdict()
    mp = multiprocessing.Pool(int(args.threads))

    for species_id, aln_stats in mp.starmap(species_pileup, species_ids):
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
        species_alnstats[species_id] = sp_stats

    return species_alnstats


def write_snps_summary(species_alnstats, outfile):
    """ Get summary of mapping statistics """
    header = ['species_id', 'genome_length', 'covered_bases', 'fraction_covered', 'mean_coverage', 'aligned_reads', 'mapped_reads']
    with OutputStream(outfile) as file:
        file.write('\t'.join(header) + '\n')
        for species_id, species_aln in species_alnstats.items():
            ## to make sure the dict key orders are preserved
            file.write('\t'.join([species_id, map(str, species_aln.values())]))
            ## moved the DECIMALS to the calculation of the values


def midas_run_snps(args):

    tempdir = f"{args.outdir}/snps/temp_sc{args.species_cov}"

    if args.debug and os.path.exists(tempdir):
        tsprint(f"INFO:  Reusing existing temp data in {tempdir} according to --debug flag.")
    else:
        command(f"rm -rf {tempdir}")
        command(f"mkdir -p {tempdir}")

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
        contigs, db_stats = initialize_contigs(contigs_files)

        # print out database stats
        tsprint(f"CONTIG DB stats - total genomes: {db_stats['species']}")
        tsprint(f"CONTIG DB stats - total contigs: {db_stats['total_seqs']}")
        tsprint(f"CONTIG DB stats - total base-pairs: {db_stats['total_length']}")

        # Use Bowtie2 to map reads to a representative genomes
        build_repgenome_db(tempdir, contigs_files)
        repgenome_align(args, tempdir)

        # Use mpileup to identify SNPs
        samtools_index(tempdir, args)
        print("hi")
        species_alnstats = pysam_pileup(args, species_profile.keys(), contigs)
        write_snps_summary(species_alnstats, f"{args.outdir}/snps/summary.txt")

    except:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error. Specify --debug flag to keep.")


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_run_snps(args)
