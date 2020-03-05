import json
from collections import defaultdict
import multiprocessing
from math import ceil

import numpy as np
from pysam import AlignmentFile  # pylint: disable=no-name-in-module
import Bio.SeqIO

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, num_physical_cores, InputStream, OutputStream, multiprocessing_map, download_reference, split, command
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
    subparser.add_argument('--species_profile_path',
                           dest='species_profile_path',
                           type=str,
                           metavar="CHAR",
                           help=f"Species profile path for the prebuild bowtie2 index")

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


def keep_read(aln):
    global global_args
    args = global_args

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

    return True


def scan_contigs(contig_file, species_id):
    contigs = {}
    with InputStream(contig_file) as file:
        for rec in Bio.SeqIO.parse(file, 'fasta'):
            contigs[rec.id] = {
                "species_id": species_id,
                "contig_id": rec.id,
                "contig_len": int(len(rec.seq)),
                "contig_seq": str(rec.seq),
            }
    return contigs


def cat_files(sliced_files, one_file, slice_num=20):
    for temp_files in split(sliced_files, slice_num):
        command("cat " + " ".join(temp_files) + f" >> {one_file}")


def merge_sliced_pileup_for_species(species_id):

    global species_sliced_snps_path
    global semaphore_for_species
    global slice_counts
    global global_args

    sliced_files = species_sliced_snps_path[species_id][:-1]
    merged_file = species_sliced_snps_path[species_id][-1]

    with OutputStream(merged_file) as stream:
        stream.write('\t'.join(snps_pileup_schema.keys()) + '\n')

    cat_files(sliced_files, merged_file, 20)

    for _ in range(slice_counts[species_id]):
        semaphore_for_species[species_id].release() # no deadlock

    if not global_args.debug:
        for s_file in sclies_files:
            command("rm -rf {s_file}")

    return True


def slice_pileup(packed_args):

    global semaphore_for_species
    global slice_counts
    global global_args
    args = global_args

    if packed_args[1] == -1:
        species_id = packed_args[0]

        for _ in range(slice_counts[species_id]):
            semaphore_for_species[species_id].acquire()

        flag = merge_sliced_pileup_for_species(species_id)
        assert flag == True, f"Failed to merge contigs snps files for species {species_id}"

    try:
        species_id, slice_id, contig_id, contig_start, contig_end, repgenome_bamfile, headerless_sliced_path, contig = packed_args
        #species_id, contig_id, repgenome_bamfile, contig, headerless_contigs_pileup_path = packed_args

        zero_rows_allowed = not args.sparse
        with AlignmentFile(repgenome_bamfile) as bamfile:
            counts = bamfile.count_coverage(
                    contig_id,
                    #start=0,
                    #end=contig["contig_len"],
                    start=contig_start,
                    end=contig_end+1,
                    quality_threshold=args.aln_baseq, # min_quality_threshold a base has to reach to be counted.
                    read_callback=keep_read) # select a call-back to ignore reads when counting

            aligned_reads = bamfile.count(
                    contig_id,
                    #start=0,
                    #end=contig["contig_len"],
                    start=contig_start,
                    end=contig_end+1)

            mapped_reads = bamfile.count(
                    contig_id,
                    #start=0,
                    #end=contig["contig_len"],
                    start=contig_start,
                    end=contig_end+1,
                    read_callback=keep_read)

        aln_stats = {
                "species_id": species_id,
                "contig_id": contig_id,
                "slice_length": contig_end - contig_start,
                #"contig_length": contig["contig_len"],
                "aligned_reads": aligned_reads,
                "mapped_reads": mapped_reads,
                "contig_total_depth": 0,
                "contig_covered_bases":0
            }

        records = []
        for ref_pos in range(contig_start, contig_end):
            print(ref_pos, contig_start, contig_end)
            ref_allele = contig["contig_seq"][ref_pos]
            depth = sum([counts[nt][ref_pos-1] for nt in range(4)])
            count_a = counts[0][ref_pos]
            count_c = counts[1][ref_pos]
            count_g = counts[2][ref_pos]
            count_t = counts[3][ref_pos]
            row = (contig_id, ref_pos+1, ref_allele, depth, count_a, count_c, count_g, count_t)

            aln_stats["contig_total_depth"] += depth
            if depth > 0:
                aln_stats["contig_covered_bases"] += 1
            if depth > 0 or zero_rows_allowed:
                records.append(row)

        with OutputStream(headerless_contigs_pileup_path) as stream:
            for row in records:
                stream.write("\t".join(map(format_data, row)) + "\n")
        return aln_stats
    finally:
        semaphore_for_species[species_id].release() # no deadlock


def species_pileup(species_ids, contigs_files, repgenome_bamfile):
    # TODO: need to provide a way to generate the temp file for each contig's SNP pileup.
    # snps_pileup_path = sample.get_target_layout("snps_pileup", species_id)
    """ Read in contigs information for one species_id """

    global global_args
    args = global_args
    global sample

    global slice_counts
    global semaphore_for_species
    global species_sliced_snps_path

    slice_counts = dict()
    semaphore_for_species = dict()
    species_sliced_snps_path = defaultdict(list)

    argument_list = []
    slice_size = 1000
    for species_index, species_id in enumerate(species_ids):
        # For each species
        contigs = scan_contigs(contigs_files[species_index], species_id)

        slice_id = 0
        for contig_id in sorted(list(contigs.keys())): # why need to sort?
            contig = contigs[contig_id]
            contig_length = contig["contig_len"]

            if contig_length <= slice_size:
                headerless_sliced_path = sample.get_target_layout("contigs_pileup", species_id, slice_id)
                slice_args = (species_id, slice_id, contig_id, 0, contig_length, repgenome_bamfile, headerless_sliced_path, contig)

                argument_list.append(slice_args)
                species_sliced_snps_path[species_id].append(headerless_sliced_path)
                slice_id += 1
            else:
                slice_num = ceil(contig_length/slice_size) - 1
                for ni, ci in enumerate(range(0, contig_length, slice_size)):
                    headerless_sliced_path = sample.get_target_layout("contigs_pileup", species_id, slice_id)

                    if ni == slice_num:
                        slice_args = (species_id, slice_id, contig_id, ci, contig_length, repgenome_bamfile, headerless_sliced_path, contig)
                    else:
                        slice_args = (species_id, slice_id, contig_id, ci, ci+slice_size, repgenome_bamfile, headerless_sliced_path, contig)

                    argument_list.append(slice_args)
                    species_sliced_snps_path[species_id].append(headerless_sliced_path)
                    slice_id += 1

        # Submit the merge jobs
        argument_list.append((species_id, -1))
        species_sliced_snps_path[species_id].append(sample.get_target_layout("snps_pileup", species_id))

        # Create a semaphore with contig_counter of elements
        semaphore_for_species[species_id] = multiprocessing.Semaphore(slice_id)
        print(slice_id)
        for _ in range(slice_id):
            semaphore_for_species[species_id].acquire()
        slice_counts[species_id] = slice_id

    contigs_pileup_summary = multiprocessing_map(slice_pileup, argument_list, num_procs=num_physical_cores)

    return contigs_pileup_summary


def compute_species_pileup_summary(contigs_pileup_summary):

    per_species_pileup_stats = {
        "genome_length": 0,
        "covered_bases": 0,
        "total_depth": 0,
        "aligned_reads":0,
        "mapped_reads":0,
        "fraction_covered": 0.0,
        "mean_coverage": 0.0
        }

    species_pileup_summary = defaultdict(dict)
    current_species_id = None

    for record in contigs_pileup_summary:

        species_id = record["species_id"]
        if species_id not in species_pileup_summary:
            species_pileup_summary[species_id] = per_species_pileup_stats

        perspecies_pileup = species_pileup_summary.get(species_id)
        perspecies_pileup["genome_length"] +=  record["species_id"]["contig_length"]
        perspecies_pileup["total_depth"] += record["species_id"]["contig_total_depth"]
        perspecies_pileup["covered_bases"] += record["species_id"]["contig_covered_bases"]
        perspecies_pileup["aligned_reads"] += record["species_id"]["aligned_reads"]
        perspecies_pileup["mapped_reads"] += record["species_id"]["mapped_reads"]

        if current_species_id is not None and current_species_id != species_id:
            current_species_pileup = species_pileup_summary.get(current_species_id)

            if current_species_pileup["genome_length"] > 0:
                current_species_pileup["fraction_covered"] = current_species_pileup["covered_bases"] / current_species_pileup["covered_bases"]
            if current_species_pileup["covered_bases"] > 0:
                current_species_pileup["mean_coverage"] = current_species_pileup["total_depth"] / current_species_pileup["covered_bases"]

            current_species_id = species_id

    if perspecies_pileup["genome_length"] > 0:
        perspecies_pileup["fraction_covered"] = perspecies_pileup["covered_bases"] / perspecies_pileup["covered_bases"]
    if perspecies_pileup["covered_bases"] > 0:
        perspecies_pileup["mean_coverage"] = perspecies_pileup["total_depth"] / perspecies_pileup["covered_bases"]

    return species_pileup_summary


def write_species_pileup_summary(species_pileup_summary, outfile):
    """ Get summary of mapping statistics """
    with OutputStream(outfile) as file:
        file.write('\t'.join(snps_profile_schema.keys()) + '\n')
        for species_id, species_summary in species_pileup_summary.items():
            record = [species_id] + list(species_summary.values())
            file.write("\t".join(map(format_data, record)) + "\n")


def midas_run_snps(args):

    global sample
    sample = Sample(args.sample_name, args.midas_outdir, "snps")
    sample.create_output_dir(args.debug)

    global global_args
    global_args = args

    try:
        if args.bt2_db_indexes:
            assert bowtie2_index_exists(bt2_db_dir, bt2_db_name) and os.path.exists(args.species_profile), f"Check the path bt2_db_dir and exists of species_profile_path"
            tsprint("Prebuild repgenomes bowtie2 index and species_profile exit. Use them")

            bt2_db_dir = os.path.dirname(args.bt2_db_indexes)
            bt2_db_name = os.path.basename(args.bt2_db_indexes)

            species_profile = {}
            with InputStream() as stream:
                for species_id, coverage in select_from_tsv(stream, ["species_id", "coverage"]):
                    species_profile[species_id] = coverage

            sample.create_species_subdir(species_profile.keys(), args.debug, "dbs")

        else:
            bt2_db_dir = sample.get_target_layout("dbsdir")
            bt2_db_name = "repgenomes"
            bt2_db_temp_dir = sample.get_target_layout("dbs_tempdir")

            if args.species_list:
                species_profile = sample.select_species(args.genome_coverage, args.species_list)
            else:
                species_profile = sample.select_species(args.genome_coverage)

            local_toc = download_reference(outputs.genomes, bt2_db_dir)
            db = UHGG(local_toc)

            # Download repgenome_id.fna for every species in the restricted species profile
            sample.create_species_subdir(species_profile.keys(), args.debug, "dbs")
            contigs_files = db.fetch_contigs(species_profile.keys(), bt2_db_temp_dir)

            build_bowtie2_db(bt2_db_dir, bt2_db_name, contigs_files)

        # Use Bowtie2 to map reads to a representative genomes
        repgenome_bamfile = sample.get_target_layout("snps_repgenomes_bam")
        bowtie2_align(bt2_db_dir, bt2_db_name, repgenome_bamfile, args)
        samtools_index(repgenome_bamfile, args.debug)


        # Create species subdir at one place
        species_ids_of_interest = species_profile.keys()
        sample.create_species_subdir(species_ids_of_interest, args.debug, "temp")
        # Use mpileup to identify SNPs
        contigs_pileup_summary = species_pileup(species_ids_of_interest, contigs_files, repgenome_bamfile)
        species_pileup_summary = compute_species_pileup_summary(contigs_pileup_summary)
        write_species_pileup_summary(species_pileup_summary, filepath)

        # TODO: for the provided bowtie2 database index, we don't have the contigs files
        # One way is to move the scan_contigs_file_stats into build database and when provicded with bowtie2 index,
        # also need to provide something like genome_stats

    except:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error. Specify --debug flag to keep.")
            sample.remove_output_dir()
        raise


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_run_snps(args)
