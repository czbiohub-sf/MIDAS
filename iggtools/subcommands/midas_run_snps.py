import json
from collections import defaultdict
import multiprocessing

import numpy as np
from pysam import AlignmentFile  # pylint: disable=no-name-in-module
import Bio.SeqIO

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, num_physical_cores, InputStream, OutputStream, multithreading_hashmap, download_reference
from iggtools.params import outputs
from iggtools.models.uhgg import UHGG
from iggtools.common.bowtie2 import build_bowtie2_db, bowtie2_align, samtools_index, bowtie2_index_exists
from iggtools.params.schemas import snps_profile_schema, snps_pileup_schema, snps_info_schema, DECIMALS, format_data
from iggtools.models.sample import Sample

DEFAULT_ALN_COV = 0.75
DEFAULT_GENOME_COVERAGE = 3.0
DEFAULT_ALN_MAPID = 94.0
DEFAULT_ALN_MAPQ = 20
DEFAULT_ALN_READQ = 20
DEFAULT_ALN_BASEQ = 30
DEFAULT_ALN_TRIM = 0


def register_args(main_func):
    subparser = add_subcommand('midas_run_snps', main_func, help='single-nucleotide-polymorphism prediction')
    subparser.add_argument('midas_outdir',
                           type=str,
                           help="""Path to directory to store results.  Name should correspond to unique sample identifier.""")
    subparser.add_argument('--sample_name',
                           dest='sample_name',
                           required=True,
                           help="Unique sample identifier")
    subparser.add_argument('-1',
                           dest='r1',
                           required=True,
                           help="FASTA/FASTQ file containing 1st mate if using paired-end reads.  Otherwise FASTA/FASTQ containing unpaired reads.")
    subparser.add_argument('-2',
                           dest='r2',
                           help="FASTA/FASTQ file containing 2nd mate if using paired-end reads.")

    subparser.add_argument('--sparse',
                           action='store_true',
                           default=False,
                           help=f"Omit zero rows from output.")

    subparser.add_argument('--genome_coverage',
                           type=float,
                           dest='genome_coverage',
                           metavar='FLOAT',
                           default=DEFAULT_GENOME_COVERAGE,
                           help=f"Include species with >X coverage ({DEFAULT_GENOME_COVERAGE})")
    # first add species_list
    subparser.add_argument('--species_list',
                           dest='species_list',
                           type=str,
                           metavar="CHAR",
                           help=f"Comma separated list of species ids")
    subparser.add_argument('--bt2_db_indexes',
                           dest='bt2_db_indexes',
                           type=str,
                           metavar="CHAR",
                           help=f"Prebuilt bowtie2 database indexes")
    # The question is: where should I put the built database ?
    # I think it's should be put at the midas_outdir <=

    # then think about midas_merge_species to get on repgenomes database
    subparser.add_argument('--pool_species',
                           dest='pool_species',
                           default=False,
                           type=str,
                           metavar="CHAR",
                           help=f"Pool species present in all samples in the given sample_list.")

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
    subparser.add_argument('--aln_sort',
                           action='store_true',
                           default=True,
                           help=f"Sort BAM file.")

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
    subparser.add_argument('--max_reads',
                           dest='max_reads',
                           type=int,
                           metavar="INT",
                           help=f"Number of reads to use from input file(s).  (All)")
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
    return main_func


def keep_read_worker(aln, pileup_stats):
    global global_args
    args = global_args

    pileup_stats['aligned_reads'] += 1

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

    pileup_stats['mapped_reads'] += 1
    return True


def species_pileup_old(species_id, repgenome_bamfile, snps_pileup_path, contig_file, contigs_db_stats):
    """ Read in contigs information for one species_id """

    global global_args
    args = global_args

    contigs_db_stats['species_counts'] += 1 # not being updated and passed as expected
    contigs = {}
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
    pileup_stats = {
        "genome_length": 0,
        "total_depth": 0,
        "covered_bases": 0,
        "aligned_reads":0,
        "mapped_reads":0,
        }

    def keep_read(x):
        return keep_read_worker(x, pileup_stats)

    with OutputStream(snps_pileup_path) as file:
        file.write('\t'.join(snps_pileup_schema.keys()) + '\n')
        zero_rows_allowed = not args.sparse

        # Loop over alignment for current species's contigs
        with AlignmentFile(repgenome_bamfile) as bamfile:
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

                    pileup_stats['genome_length'] += 1
                    pileup_stats['total_depth'] += depth
                    if depth > 0:
                        pileup_stats['covered_bases'] += 1

    tsprint(json.dumps({species_id: pileup_stats}, indent=4))
    return (species_id, {k: str(v) for k, v in pileup_stats.items()})


def species_pileup(species_id, repgenome_bamfile, snps_pileup_path, contig_file, contigs_db_stats):
    """ Read in contigs information for one species_id """

    global global_args
    args = global_args

    contigs_db_stats['species_counts'] += 1 # not being updated and passed as expected
    # CONTIGS_DB_STATS didn't work as expected, so I remove codes related to that for now. Add me back later on after figuring out the parallel.

    contigs = {}
    with InputStream(contig_file) as file:
        for rec in Bio.SeqIO.parse(file, 'fasta'):
            contigs[rec.id] = {
                "species_id": species_id,
                "contig_len": int(len(rec.seq)),
                "contig_seq": str(rec.seq),
            }


    genome_length = sum(c.contig_len for c in contigs)
    print("genome_length", genome_length)

    # Summary statistics
    pileup_stats = {
        "genome_length": 0, # we already know this number, isn't it?
        "total_depth": 0,
        "covered_bases": 0,
        "aligned_reads":0,
        "mapped_reads":0,
        }

    def keep_read(x):
        return keep_read_worker(x, pileup_stats)

    with OutputStream(snps_pileup_path) as file:
        file.write('\t'.join(snps_pileup_schema.keys()) + '\n')
        zero_rows_allowed = not args.sparse

        # at least we need to paralize of the unit of contigs
        # Loop over alignment for current species's contigs
        #with AlignmentFile(repgenome_bamfile) as bamfile:
        for contig_id in sorted(list(contigs.keys())): # why need to sort?

            contig = contigs[contig_id]

            pileup_stats, contigs, repgenome_bamfile, contig_id  =  packged_args
            contig = contigs[contig_id]

            contig_pileup = defaultdict()

            # through this value we can easily know the ref_allele without need to loop over the BAM file
            #(contig["contig_seq"][ref_pos] for ref_pos in range(0, contig["contig_len"]))

def process_contig(packged_args):
    pileup_stats, contigs, repgenome_bamfile, contig_id = packged_args

    global global_args
    args = global_args
    results = []
    with AlignmentFile(repgenome_bamfile) as bamfile:
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
        record = (contig_id, ref_pos + 1, ref_allele, depth, count_a, count_c, count_g, count_t)

        # Do I write separately or do I write to the same file by multiprocessingself
        # Write to the same file while compute by different processes seems definitely a bad idea.
        if depth > 0 or zero_rows_allowed:
            #file.write("\t".join(map(format_data, record)) + "\n")
            results.append(record)

        pileup_stats['genome_length'] += 1 # why do we do this here??
        pileup_stats['total_depth'] += depth
        if depth > 0:
            pileup_stats['covered_bases'] += 1

    tsprint(json.dumps({species_id: pileup_stats}, indent=4))
    return (species_id, {k: str(v) for k, v in pileup_stats.items()})



def pysam_pileup_old(species_ids, contigs_files):
    """ Counting alleles and run pileups per species in parallel """
    global global_sample
    sample = global_sample

    # Update alignment stats for species
    species_pileup_stats = defaultdict()
    contigs_db_stats = {'species_counts':0, 'total_seqs':0, 'total_length':0}

    argument_list = []
    for species_index, species_id in enumerate(species_ids):
        my_args = (species_id, sample.get_target_layout("snps_repgenomes_bam"), sample.get_target_layout("snps_pileup", species_id), contigs_files[species_index], contigs_db_stats)
        argument_list.append(my_args)

    mp = multiprocessing.Pool(num_physical_cores)
    for species_id, pileup_stats in mp.starmap(species_pileup, argument_list):
        species_stats = {
            "genome_length": int(pileup_stats['genome_length']),
            "covered_bases": int(pileup_stats['covered_bases']),
            "total_depth": int(pileup_stats['total_depth']),
            "aligned_reads": int(pileup_stats['aligned_reads']),
            "mapped_reads": int(pileup_stats['mapped_reads']),
            "fraction_covered": 0.0,
            "mean_coverage": 0.0,
        }
        if species_stats["genome_length"] > 0:
            species_stats["fraction_covered"] = species_stats["covered_bases"] / species_stats["genome_length"]
        if species_stats["covered_bases"] > 0:
            species_stats["mean_coverage"] = species_stats["total_depth"] / species_stats["covered_bases"]

        species_pileup_stats[species_id] = species_stats

    tsprint(f"contigs_db_stats - total genomes: {contigs_db_stats['species_counts']}")
    tsprint(f"contigs_db_stats - total contigs: {contigs_db_stats['total_seqs']}")
    tsprint(f"contigs_db_stats - total base-pairs: {contigs_db_stats['total_length']}")

    return species_pileup_stats




def write_snps_summary(species_pileup_stats, outfile):
    """ Get summary of mapping statistics """
    with OutputStream(outfile) as file:
        file.write('\t'.join(snps_profile_schema.keys()) + '\n')
        for species_id, species_aln in species_pileup_stats.items():
            values = list(species_aln.values())
            values.insert(0, species_id)
            file.write('\t'.join(map(str, values)) + '\n')


def midas_run_snps(args):

    sample = Sample(args.sample_name, args.midas_outdir, "snps")
    sample.create_output_dir(args.debug)

    global global_sample
    global_sample = sample
    global global_args
    global_args = args

    try:
        if args.bt2_db_indexes:
            if bowtie2_index_exists(os.path.dirname(args.bt2_db_indexes), os.path.basename(bt2_db_indexes)):
                print("good the pre built bt2 index exist")
                bt2_db_dir = os.path.dirname(args.bt2_db_indexes)
                bt2_db_name = os.path.basename(bt2_db_indexes)
            else:
                print("Error: good the pre built repgrenomes bt2 index exist")
                exit(0)
        else:
            bt2_db_dir = sample.get_target_layout("dbsdir")
            bt2_db_temp_dir = sample.get_target_layout("dbs_tempdir")
            bt2_db_name = "repgenomes"

            if args.species_list:
                species_profile = sample.select_species(args.genome_coverage, args.species_list)
            else:
                species_profile = sample.select_species(args.genome_coverage)

            local_toc = download_reference(outputs.genomes, bt2_db_dir)
            db = UHGG(local_toc)

            # Download repgenome_id.fna for every species in the restricted species profile
            sample.create_species_subdir(species_profile.keys())
            contigs_files = db.fetch_contigs(species_profile.keys(), bt2_db_temp_dir)

            build_bowtie2_db(bt2_db_dir, bt2_db_name, contigs_files)

        # Use Bowtie2 to map reads to a representative genomes
        bowtie2_align(bt2_db_dir, bt2_db_name, args)

        # Use mpileup to identify SNPs
        samtools_index(args, bt2_db_dir, bt2_db_name)
        species_pileup_stats = pysam_pileup(species_profile.keys(), contigs_files)
        # TODO: for the provided bowtie2 database index, we don't have the contigs files
        # One way is to move the scan_contigs_file_stats into build database and when provicded with bowtie2 index,
        # also need to provide something like genome_stats
        write_snps_summary(species_pileup_stats, sample.get_target_layout("snps_summary"))

    except:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error. Specify --debug flag to keep.")
            sample.remove_output_dir()


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_run_snps(args)
