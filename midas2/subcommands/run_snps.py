#!/usr/bin/env python3
import json
import os
import multiprocessing
from math import floor
from operator import itemgetter
from collections import defaultdict

import numpy as np
from pysam import AlignmentFile  # pylint: disable=no-name-in-module

from midas2.common.argparser import add_subcommand
from midas2.common.utils import tsprint, InputStream, OutputStream, multiprocessing_map, command, cat_files, select_from_tsv, multithreading_map
from midas2.common.bowtie2 import build_bowtie2_db, bowtie2_align, samtools_sort, samtools_index, bowtie2_index_exists, _keep_read
from midas2.params.schemas import snps_profile_schema, snps_pileup_schema, format_data, snps_chunk_summary_schema, snps_pileup_basic_schema
from midas2.common.snvs import call_alleles, reference_overlap, update_overlap, mismatches_within_overlaps, query_overlap_qualities
from midas2.common.utilities import scan_fasta
from midas2.models.midasdb import MIDAS_DB
from midas2.models.sample import Sample
from midas2.models.species import Species, parse_species
from midas2.params.inputs import MIDASDB_NAMES


DEFAULT_MARKER_DEPTH = 5.0
DEFAULT_MARKER_MEDIAN_DEPTH = 2

DEFAULT_ALN_MAPID = 94.0
DEFAULT_ALN_MAPQ = 10
DEFAULT_ALN_COV = 0.75
DEFAULT_ALN_READQ = 20
DEFAULT_ALN_BASEQ = 30
DEFAULT_ALN_TRIM = 0

DEFAULT_CHUNK_SIZE = 1000000
DEFAULT_MAX_FRAGLEN = 5000
DEFAULT_NUM_CORES = 8

DEFAULT_SITE_DEPTH = 2
DEFAULT_SNP_MAF = 0.1


def register_args(main_func):
    subparser = add_subcommand('run_snps', main_func, help='Predict single-nucleotide-polymorphism')

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

    # Prebuilt/predownload
    subparser.add_argument('--prebuilt_bowtie2_indexes',
                           dest='prebuilt_bowtie2_indexes',
                           type=str,
                           metavar="CHAR",
                           help=f"Prebuilt bowtie2 database indexes")
    subparser.add_argument('--prebuilt_bowtie2_species',
                           dest='prebuilt_bowtie2_species',
                           type=str,
                           metavar="CHAR",
                           help=f"List of species used for building the prebuild bowtie2 indexes.")
    subparser.add_argument('--midasdb_name',
                           dest='midasdb_name',
                           type=str,
                           default="uhgg",
                           choices=MIDASDB_NAMES,
                           help=f"MIDAS Database name.")
    subparser.add_argument('--midasdb_dir',
                           dest='midasdb_dir',
                           type=str,
                           default="midasdb",
                           help=f"Local MIDAS Database path mirroing S3.")

    # Species related
    subparser.add_argument('--species_list',
                           dest='species_list',
                           type=str,
                           metavar="CHAR",
                           help=f"Comma separated list of species ids")
    subparser.add_argument('--select_by',
                           dest='select_by',
                           type=str,
                           default="median_marker_coverage",
                           #choices=['median_marker_coverage', 'marker_coverage', 'unique_fraction_covered', "marker_relative_abundance"],
                           help=f"Comma separated columns from species_profile to filter species.")
    subparser.add_argument('--select_threshold',
                           dest='select_threshold',
                           type=str,
                           metavar="CHAR",
                           default=str(DEFAULT_MARKER_MEDIAN_DEPTH),
                           help=f"Comman separated correponsding cutoff for filtering species (> XX) ({DEFAULT_MARKER_MEDIAN_DEPTH}, )")

    # Alignment flags (Bowtie2 or postprocessing)
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
                           help='Global/local read alignment (local, corresponds to the bowtie2 --local; Default global corresponds to the bowtie2 --end-to-end).')
    subparser.add_argument('--aln_interleaved',
                           action='store_true',
                           default=False,
                           help='FASTA/FASTQ file in -1 are paired and contain forward AND reverse reads')
    subparser.add_argument('--fragment_length',
                           type=float,
                           dest='fragment_length',
                           metavar='FLOAT',
                           default=DEFAULT_MAX_FRAGLEN,
                           help=f"Maximum fragment length for paired reads ({DEFAULT_MAX_FRAGLEN}).")

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

    subparser.add_argument('--paired_only',
                           action='store_true',
                           default=False,
                           help=f"Only recruit properly paired reads for pileup.")
    subparser.add_argument('--advanced',
                           action='store_true',
                           default=False,
                           help=f"Report majore/minor allele for each genomic sites.")
    subparser.add_argument('--analysis_ready',
                           action='store_true',
                           default=False,
                           help=f"Report majore/minor allele for each genomic sites.")

    # File related
    subparser.add_argument('--chunk_size',
                           dest='chunk_size',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_CHUNK_SIZE,
                           help=f"Number of genomic sites for the temporary chunk file  ({DEFAULT_CHUNK_SIZE})")
    subparser.add_argument('--max_reads',
                           dest='max_reads',
                           type=int,
                           metavar="INT",
                           help=f"Number of reads to use from input file(s).  (All)")
    subparser.add_argument('--num_cores',
                           dest='num_cores',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_NUM_CORES,
                           help=f"Number of physical cores to use ({DEFAULT_NUM_CORES})")

    subparser.add_argument('--site_depth',
                           dest='site_depth',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_SITE_DEPTH,
                           help=f"Minimum number of reads mapped to genomic site ({DEFAULT_SITE_DEPTH})")
    subparser.add_argument('--snp_maf',
                           dest='snp_maf',
                           type=float,
                           metavar="FLOAT",
                           default=DEFAULT_SNP_MAF,
                           help=f"Minimum single sample minor-allele_frequency to call an allele present ({DEFAULT_SNP_MAF}), Values > 0.0 and < 0.5 are accepted.")

    return main_func


def keep_read(aln):
    global global_args
    args = global_args
    if not args.paired_only:
        return _keep_read(aln, args.aln_mapid, args.aln_readq, args.aln_mapq, args.aln_cov)
    return True


def in_place(species_counts):
    # Carry contigs in global variable
    return species_counts < 50


def design_chunks_per_species(args):
    sp, midas_db, chunk_size, carry = args
    if carry:
        sp.get_repgenome(midas_db)
    return sp.compute_snps_chunks(midas_db, chunk_size, "run")


def design_chunks(species_ids_of_interest, midas_db, chunk_size):
    """ Chunks of continuous genomics sites, indexed by species_id, chunk_id """

    global semaphore_for_species
    global dict_of_species
    global dict_of_site_chunks

    # Read-only global variables
    semaphore_for_species = dict()
    dict_of_species = {species_id: Species(species_id) for species_id in species_ids_of_interest}

    # Design chunks structure per species
    num_cores = min(midas_db.num_cores, 16)
    carry = in_place(len(dict_of_species))
    all_site_chunks = multithreading_map(design_chunks_per_species, [(sp, midas_db, chunk_size, carry) for sp in dict_of_species.values()], num_cores) #<---

    dict_of_site_chunks = defaultdict(dict)
    for spidx, species_id in enumerate(species_ids_of_interest):
        dict_of_site_chunks[species_id] = all_site_chunks[spidx]

    # Sort species by the max_contig_length or num_of_snps_chunks
    sorted_tuples_of_species = sorted(((sp.id, sp.num_of_snps_chunks) for sp in dict_of_species.values()), key=itemgetter(1), reverse=True)

    arguments_list = []
    for species_id, _ in sorted_tuples_of_species:
        sp = dict_of_species[species_id]

        num_of_snps_chunks = sp.num_of_snps_chunks
        for chunk_id in range(0, num_of_snps_chunks):
            arguments_list.append((species_id, chunk_id))

        # Create a semaphore with number_of_chunks of elements
        semaphore_for_species[species_id] = multiprocessing.Semaphore(num_of_snps_chunks)
        for _ in range(num_of_snps_chunks):
            semaphore_for_species[species_id].acquire()

    tsprint("================= Total number of compute chunks: " + str(len(arguments_list)))

    for species_id in dict_of_species.keys():
        arguments_list.append((species_id, -1))

    return arguments_list


def filter_bam_by_single_read(species_id, repbamfile, filtered_bamfile):
    """ Filter given BAM file with propely paired reads for given species and write to file """

    global global_args
    global dict_of_species
    global sample

    tsprint(f"  MIDAS::filter_bam_by_single_read::{species_id}-0::start filter_bam_by_single_read")

    # List of contigs for given species
    sp = dict_of_species[species_id]

    if in_place(len(dict_of_species)):
        list_of_contig_ids = list(dict_of_species[species_id].contigs.keys())
    else:
        list_of_contig_ids = sp.fetch_contigs_ids()

    # Cache *properly* aligned reads-pair
    filtered_alns_dict = defaultdict(dict)
    reads_stats = {
        "aligned_reads": dict.fromkeys(list_of_contig_ids, 0),
        "mapped_reads": dict.fromkeys(list_of_contig_ids, 0)
    }

    with AlignmentFile(repbamfile) as infile:
        for contig_id in list_of_contig_ids:
            # To avoid overcount boudary reads, we compute reads stats per contig.
            aligned_reads = 0
            mapped_reads = 0
            for aln in infile.fetch(contig_id):
                aligned_reads += 1
                if global_args.analysis_ready or keep_read(aln):
                    mapped_reads += 1
                    read = "1" if aln.is_read1 else "2"
                    filtered_alns_dict[f"{aln.query_name}_{read}"] = aln
            reads_stats["aligned_reads"][contig_id] = aligned_reads
            reads_stats["mapped_reads"][contig_id] = mapped_reads

    # Write filtered alignments to file
    template_bam = AlignmentFile(repbamfile, "rb")
    filtered_bam = AlignmentFile(filtered_bamfile, "wb", template=template_bam)
    for aln in filtered_alns_dict.values():
        filtered_bam.write(aln)
    filtered_bam.close()

    tsprint(f"  MIDAS::filter_bam_by_single_read::{species_id}-0::finish filter_bam_by_single_read")
    return reads_stats


def filter_bam_by_proper_pair(species_id, repbamfile, filtered_bamfile):
    """ Filter given BAM file with propely paired reads for given species and write to file """

    global global_args
    global dict_of_species
    global sample

    tsprint(f"  MIDAS::filter_bam_by_proper_pair::{species_id}-0::start filter_bam_by_proper_pair")

    # List of contigs for given species
    sp = dict_of_species[species_id]

    if in_place(len(dict_of_species)):
        list_of_contig_ids = list(dict_of_species[species_id].contigs.keys())
    else:
        list_of_contig_ids = sp.fetch_contigs_ids()

    # Cache *properly* aligned reads-pair
    filtered_alns_dict = defaultdict(dict)
    reads_stats = {
        "aligned_reads": dict.fromkeys(list_of_contig_ids, 0),
        "mapped_reads": dict.fromkeys(list_of_contig_ids, 0)
    }

    with AlignmentFile(repbamfile) as infile:
        for contig_id in list_of_contig_ids:
            # To avoid boundary cliff, we need to read in the whole contig
            aligned_reads = 0
            alns_dict = defaultdict(dict) # cache the reads
            for aln in infile.fetch(contig_id):
                aligned_reads += 1
                if aln.is_secondary:
                    continue
                if not aln.is_proper_pair:
                    continue
                if aln.is_reverse:
                    alns_dict[aln.query_name]["rev"] = aln
                else:
                    alns_dict[aln.query_name]["fwd"] = aln
            reads_stats["aligned_reads"][contig_id] = aligned_reads

            mapped_reads = 0
            for query_name, alns in alns_dict.items():
                # Ignore orphan reads
                if len(alns) != 2:
                    continue

                # Apply filters to paired-reads
                # Common features
                readq = np.mean(alns["fwd"].query_qualities + alns["rev"].query_qualities)
                mapq = max(alns["fwd"].mapping_quality, alns["rev"].mapping_quality)
                if readq < global_args.aln_readq:
                    continue
                if mapq < global_args.aln_mapq:
                    continue

                # Template length: number of bases from the left most mapped base to the rightmost mapped base on the reference
                fragment_length = abs(alns["fwd"].template_length)
                if fragment_length >= global_args.fragment_length:
                    continue

                # I think the alignment coverage should not be affected by overlap.
                # However, we should double check whether gaps counted as aligned ..
                align_len = alns["fwd"].query_alignment_length + alns["rev"].query_alignment_length
                query_len = alns["fwd"].query_length + alns["rev"].query_length
                alncov = align_len / float(query_len)
                if alncov < global_args.aln_cov:
                    continue

                # For the compute of sequence identity, we need to specially consider paired-reads overlap
                # Compute the length of the overlapping region along the reference
                reads_overlap = reference_overlap((alns["fwd"].reference_start, alns["fwd"].reference_end - 1), (alns["rev"].reference_start, alns["rev"].reference_end - 1))
                # Compute the query overlap length: substract the gaps in the aligned from the FWD reads to define the overlap boundary
                reads_overlap = update_overlap(reads_overlap, alns["fwd"])

                if reads_overlap:
                    # Keep the FWD read, split the REV reads
                    (nm_out_rev, nm_in_rev, _, _) = mismatches_within_overlaps(alns["rev"], reads_overlap, "rev")
                    #assert nm_out_rev + nm_in_rev == dict(alns["rev"].tags)['NM'], f"REV {query_name}"

                    # Keep the REV read, split the FWD reads
                    (nm_out_fwd, nm_in_fwd, ngaps_ri_fwd, _) = mismatches_within_overlaps(alns["fwd"], reads_overlap, "fwd")
                    #assert nm_out_fwd + nm_in_fwd == dict(alns["fwd"].tags)['NM'], f"FWD {query_name}"

                    # Update the overlap by substracting the number of gaps in the fwd overlap region
                    reads_overlap = reads_overlap - ngaps_ri_fwd

                    # For repeats regions, paired-end reads can be aligned with many gaps, resulting in high mismatches within the overlapping region
                    # Only keep aligned pairs indicating from the same DNA fragment
                    if abs(nm_in_fwd - nm_in_rev) > 1:
                        continue #<-----------


                    mismatches = dict(alns["fwd"].tags)['NM'] + nm_out_rev
                    align_len = alns["rev"].query_alignment_length + alns["fwd"].query_alignment_length - reads_overlap
                    mapid = 100 * (align_len - mismatches) / float(align_len)

                    # To avoid overcounting site depth for the overlapping region,
                    # "The higher quality base is used and the lower-quality base is set to BQ=0."
                    b1 = alns["fwd"].query_alignment_end - reads_overlap
                    b2 = alns["rev"].query_alignment_start + reads_overlap - 1

                    # Only use the higher quality base in the overlap region for downstream pileup
                    f = alns["fwd"].query_qualities[b1:]
                    r = alns["rev"].query_qualities[:b2+1]
                    f, r = query_overlap_qualities(f, r)
                    alns["fwd"].query_qualities[b1:] = f
                    alns["rev"].query_qualities[:b2+1] = r
                else:
                    mismatches = dict(alns["fwd"].tags)['NM'] + dict(alns["rev"].tags)['NM']
                    mapid = 100 * (align_len - mismatches) / float(align_len)

                if mapid < global_args.aln_mapid:
                    continue

                # Compute the mapped reads for the whole contig
                mapped_reads += 2
                filtered_alns_dict[query_name]["fwd"] = alns["fwd"]
                filtered_alns_dict[query_name]["rev"] = alns["rev"]
            reads_stats["mapped_reads"][contig_id] = mapped_reads

    # Write filtered alignments to file
    template_bam = AlignmentFile(repbamfile, "rb")
    filtered_bam = AlignmentFile(filtered_bamfile, "wb", template=template_bam)
    for query_name, alns in filtered_alns_dict.items():
        filtered_bam.write(alns["fwd"])
        filtered_bam.write(alns["rev"])
    filtered_bam.close()

    tsprint(f"  MIDAS::filter_bam_by_proper_pair::{species_id}-0::finish filter_bam_by_proper_pair")
    return reads_stats


def filter_bam(pargs):
    global global_args
    species_id, repgenome_bamfile, filtered_bamfile = pargs

    if global_args.paired_only:
        reads_stats = filter_bam_by_proper_pair(species_id, repgenome_bamfile, filtered_bamfile)
    else:
        reads_stats = filter_bam_by_single_read(species_id, repgenome_bamfile, filtered_bamfile)

    return reads_stats


def sort_bam(pargs):
    global global_args
    global sample

    species_id, numcores = pargs

    filtered_bamfile = sample.get_target_layout("species_bam", species_id) # <species_id> output
    sorted_bamfile = sample.get_target_layout("species_sorted_bam", species_id)

    samtools_sort(filtered_bamfile, sorted_bamfile, global_args.debug, numcores)
    samtools_index(sorted_bamfile, global_args.debug, numcores)
    command(f"rm -rf {filtered_bamfile}", quiet=True)

    return True


def process_chunk_of_sites(packed_args):
    """ Process one chunk: either pileup or merge and write results to disk """

    species_id, chunk_id = packed_args

    if chunk_id == -1:
        global semaphore_for_species
        global dict_of_species
        sp = dict_of_species[species_id]

        tsprint(f"  MIDAS::process_chunk_of_sites::{species_id}-{chunk_id}::wait merge_chunks_per_species")
        for _ in range(sp.num_of_snps_chunks):
            semaphore_for_species[species_id].acquire()
        tsprint(f"  MIDAS::process_chunk_of_sites::{species_id}-{chunk_id}::start merge_chunks_per_species")
        ret = merge_chunks_per_species(species_id)
        tsprint(f"  MIDAS::process_chunk_of_sites::{species_id}-{chunk_id}::finish merge_chunks_per_species")
        return ret

    tsprint(f"  MIDAS::process_chunk_of_sites::{species_id}-{chunk_id}::start compute_pileup_per_chunk")
    ret = compute_pileup_per_chunk(packed_args)
    tsprint(f"  MIDAS::process_chunk_of_sites::{species_id}-{chunk_id}::finish compute_pileup_per_chunk")

    return ret


def compute_pileup_per_chunk(packed_args):
    """ Pileup for one chunk, potentially contain multiple contigs """

    global semaphore_for_species
    global dict_of_species
    global sample
    global global_args
    global dict_of_site_chunks

    try:
        species_id, chunk_id = packed_args
        sp = dict_of_species[species_id]

        chunks_of_sites = dict_of_site_chunks[species_id]
        if in_place(len(dict_of_species)):
            contigs = dict_of_species[species_id].contigs
        else:
            contigs = scan_fasta(sp.contigs_fp)

        dict_of_chunk_pileup = defaultdict(dict)
        ret = []
        for pidx, pargs in enumerate(chunks_of_sites[chunk_id]):
            species_id, chunk_id, contig_id, contig_start, contig_end, count_flag = pargs[:6]

            curr_contig = contigs[contig_id]
            aln_stats, sliced_pileup = midas_pileup((species_id, chunk_id, contig_id, contig_start, contig_end, count_flag, curr_contig))
            ret.append(aln_stats)
            dict_of_chunk_pileup[pidx] = sliced_pileup

        headerless_sliced_path = sample.get_target_layout("chunk_pileup", species_id, chunk_id)
        with OutputStream(headerless_sliced_path) as stream:
            for sliced_pileup in dict_of_chunk_pileup.values():
                for row in sliced_pileup:
                    stream.write("\t".join(map(format_data, row)) + "\n")
        return ret
    finally:
        semaphore_for_species[species_id].release() # no deadlock


def midas_pileup(packed_args):
    """ Pileup for continuous positions of one contig in one chunk """

    global global_args
    global dict_of_species
    global sample

    # [contig_start, contig_end)
    species_id, chunk_id, contig_id, contig_start, contig_end, _, curr_contig = packed_args

    repgenome_bamfile = sample.get_target_layout("species_sorted_bam", species_id)

    current_chunk_size = contig_end - contig_start
    contig_seq = curr_contig["seq"]

    with AlignmentFile(repgenome_bamfile) as bamfile:
        counts = bamfile.count_coverage(contig_id, contig_start, contig_end, quality_threshold=global_args.aln_baseq)

    aligned_reads = 0
    mapped_reads = 0

    # aln_stats need to be passed from child process back to parents
    aln_stats = {
        "species_id": species_id,
        "chunk_id": chunk_id,
        "contig_id": contig_id,
        "chunk_length": current_chunk_size,
        "aligned_reads": aligned_reads,
        "mapped_reads": mapped_reads,
        "contig_total_depth": 0,
        "contig_covered_bases": 0,
    }

    sliced_pileup = list()
    for within_chunk_index in range(0, current_chunk_size):
        depth = sum([counts[nt][within_chunk_index] for nt in range(4)])
        count_a = counts[0][within_chunk_index]
        count_c = counts[1][within_chunk_index]
        count_g = counts[2][within_chunk_index]
        count_t = counts[3][within_chunk_index]

        ref_pos = within_chunk_index + contig_start
        ref_allele = contig_seq[ref_pos]
        row = [contig_id, ref_pos + 1, ref_allele, depth, count_a, count_c, count_g, count_t]

        if depth < global_args.site_depth:
            continue

        aln_stats["contig_total_depth"] += depth
        aln_stats["contig_covered_bases"] += 1

        if global_args.advanced:
            # Compuate single sample major/minor allele
            rc_ACGT = [count_a, count_c, count_g, count_t]
            tuple_of_alleles = zip(['A', 'C', 'G', 'T'], rc_ACGT)
            major_allele, minor_allele, _, allele_counts = call_alleles(tuple_of_alleles, depth, global_args.snp_maf)

            if allele_counts == 0:
                continue

            major_index = 'ACGT'.index(major_allele)
            minor_index = 'ACGT'.index(minor_allele)
            major_allelefreq = rc_ACGT[major_index] / depth
            minor_allelefreq = 0.0 if major_index == minor_index else rc_ACGT[minor_index] / depth

            row.extend([major_allele, minor_allele, major_allelefreq, minor_allelefreq, allele_counts])

        sliced_pileup.append(tuple(row)) # list of tuples_of_row_record

    assert within_chunk_index+contig_start == contig_end-1, f"compute_pileup_per_chunk::index mismatch error for {contig_id}."

    return aln_stats, sliced_pileup


def merge_chunks_per_species(species_id):
    """ merge the pileup results from chunks into one file per species """

    global global_args
    global sample
    global dict_of_species

    sp = dict_of_species[species_id]
    number_of_chunks = sp.num_of_snps_chunks

    list_of_chunks_pileup = [sample.get_target_layout("chunk_pileup", species_id, chunk_id) for chunk_id in range(0, number_of_chunks)]
    species_snps_pileup_file = sample.get_target_layout("snps_pileup", species_id)

    with OutputStream(species_snps_pileup_file) as stream:
        if global_args.advanced:
            stream.write("\t".join(snps_pileup_schema.keys()) + "\n")
        else:
            stream.write("\t".join(snps_pileup_basic_schema.keys()) + "\n")
    cat_files(list_of_chunks_pileup, species_snps_pileup_file, 20)

    if global_args.analysis_ready or not global_args.debug:
        tsprint(f"Deleting temporary sliced pileup files for {species_id}.")
        for s_file in list_of_chunks_pileup:
            command(f"rm -rf {s_file}", quiet=True)

        repgenome_bamfile = sample.get_target_layout("species_sorted_bam", species_id)
        if not global_args.analysis_ready:
            command(f"rm -rf {repgenome_bamfile}", quiet=True)
            command(f"rm -rf {repgenome_bamfile}.bai", quiet=True)

    # return a status flag
    # the path should be computable somewhere else
    return True


def compute_chunk_aln_summary(list_of_contig_aln_stats, species_ids_of_interest):
    """ Collect Compute chunk-level alignment stats from contigs' mapping summary"""
    global global_args
    global dict_of_site_chunks

    dict_of_chunk_aln_stats = dict()

    for spidx, species_id in enumerate(species_ids_of_interest):

        cc_to_ch = defaultdict(lambda: defaultdict(dict))
        chunks_of_sites = dict_of_site_chunks[species_id]

        for chunk_id, tchunks_list in chunks_of_sites.items():
            if chunk_id == -1:
                continue
            #[(species_id, chunk_id, contig_id, ci, ci+chunk_size, count_flag, 0)]
            for loc in tchunks_list:
                contig_id = loc[2]
                count_flag = loc[5]

                aligned_reads = list_of_contig_aln_stats[spidx]["aligned_reads"][contig_id]
                mapped_reads = list_of_contig_aln_stats[spidx]["mapped_reads"][contig_id]

                if contig_id not in cc_to_ch[chunk_id]:
                    cc_to_ch[chunk_id][contig_id] = {"aligned_reads": 0, "mapped_reads": 0}

                cc_to_ch[chunk_id][contig_id]["aligned_reads"] = aligned_reads if count_flag else 0
                cc_to_ch[chunk_id][contig_id]["mapped_reads"] = mapped_reads if count_flag else 0

        dict_of_chunk_aln_stats[species_id] = cc_to_ch
    return dict_of_chunk_aln_stats


def write_species_pileup_summary(chunks_pileup_summary, snps_summary_outfile, chunk_output, dict_of_chunk_aln_stats):
    """ Collect species pileup aln stats from all chunks and write to file """

    species_pileup_summary = defaultdict(dict)

    with OutputStream(chunk_output) as stream:
        stream.write("\t".join(snps_chunk_summary_schema.keys()) + "\n")

        for records in chunks_pileup_summary:
            if records is True:
                continue

            for record in records:
                species_id = record["species_id"]
                chunk_id = record["chunk_id"]
                contig_id = record["contig_id"]

                record["aligned_reads"] = dict_of_chunk_aln_stats[species_id][chunk_id][contig_id]["aligned_reads"]
                record["mapped_reads"] = dict_of_chunk_aln_stats[species_id][chunk_id][contig_id]["mapped_reads"]

                stream.write("\t".join(map(format_data, record.values())) + "\n")

                if species_id not in species_pileup_summary:
                    species_pileup_summary[species_id] = {
                        "species_id": species_id,
                        "genome_length": 0,
                        "covered_bases": 0,
                        "total_depth": 0,
                        "aligned_reads": 0,
                        "mapped_reads": 0,
                        "fraction_covered": 0.0,
                        "mean_coverage": 0.0
                        }

                curr_species_pileup = species_pileup_summary.get(species_id)
                curr_species_pileup["genome_length"] += record["chunk_length"]
                curr_species_pileup["total_depth"] += record["contig_total_depth"]
                curr_species_pileup["covered_bases"] += record["contig_covered_bases"]
                curr_species_pileup["aligned_reads"] += record["aligned_reads"]
                curr_species_pileup["mapped_reads"] += record["mapped_reads"]


    # Secondary round compute: need to loop over species to compute fraction_covered
    for species_id in species_pileup_summary.keys():
        curr_species_pileup = species_pileup_summary.get(species_id)
        if curr_species_pileup["genome_length"] > 0:
            curr_species_pileup["fraction_covered"] = curr_species_pileup["covered_bases"] / curr_species_pileup["genome_length"]
        if curr_species_pileup["covered_bases"] > 0:
            curr_species_pileup["mean_coverage"] = curr_species_pileup["total_depth"] / curr_species_pileup["covered_bases"]

    # Write to file
    with OutputStream(snps_summary_outfile) as stream:
        stream.write("\t".join(snps_profile_schema.keys()) + "\n")
        for record in species_pileup_summary.values():
            stream.write("\t".join(map(format_data, record.values())) + "\n")


def run_snps(args):

    try:
        global global_args
        global_args = args

        assert not (args.analysis_ready and args.paired_only), f"For analysis-ready BAM file, set --paired_only to False"

        global sample
        sample = Sample(args.sample_name, args.midas_outdir, "snps")
        sample.create_dirs(["outdir", "tempdir"], args.debug, quiet=True)

        species_list = parse_species(args)

        # Prepare Bowtie2 genome database path and name, prebuilt or to be built.
        if args.prebuilt_bowtie2_indexes:
            bt2db_dir = os.path.dirname(args.prebuilt_bowtie2_indexes)
            bt2db_name = os.path.basename(args.prebuilt_bowtie2_indexes)

            assert bowtie2_index_exists(bt2db_dir, bt2db_name), f"Provided {bt2db_dir}/{bt2db_name} don't exist."
            assert (args.prebuilt_bowtie2_species and os.path.exists(args.prebuilt_bowtie2_species)), f"Require list of speices used to build the provided Bowtie2 indexes." # to fetch fasta seq

            tsprint(f"Read in list of species used to build provided bowtie2 indexes {bt2db_dir}/{bt2db_name}")
            with InputStream(args.prebuilt_bowtie2_species) as stream:
                bt2_species_list = [spid[0] for spid in select_from_tsv(stream, schema={"species_id": str})]

            # Update the species list
            species_list = list(set(species_list) & set(bt2_species_list)) if species_list else bt2_species_list
        else:
            sample.create_dirs(["bt2_indexes_dir"], args.debug, quiet=True)
            bt2db_dir = sample.get_target_layout("bt2_indexes_dir")
            bt2db_name = "repgenomes"

        # Pileup species: abundant and/or listed species. We don't recommend pileup on too many empty species.
        select_thresholds = args.select_threshold.split(',')
        no_filter = len(select_thresholds) == 1 and float(select_thresholds[0]) == -1

        species_ids_of_interest = species_list if no_filter else sample.select_species(args, species_list)
        species_counts = len(species_ids_of_interest)

        sample.create_species_subdirs(species_ids_of_interest, "temp", args.debug, quiet=True)
        assert species_counts > 0, f"No (specified) species pass the marker_depth filter, please adjust the marker_depth or species_list"
        tsprint(len(species_ids_of_interest))

        # Fetch representative genome fastas for each species (multiprocessing)
        tsprint(f"MIDAS::design_chunks::start")
        num_cores_download = min(args.num_cores, species_counts)
        midas_db = MIDAS_DB(os.path.abspath(args.midasdb_dir), args.midasdb_name, num_cores_download)

        arguments_list = design_chunks(species_ids_of_interest, midas_db, args.chunk_size)
        tsprint(f"MIDAS::design_chunks::finish")

        # Build Bowtie indexes for species in the restricted species profile
        contigs_files = midas_db.fetch_files("representative_genome", species_ids_of_interest)
        tsprint(f"MIDAS::build_bowtie2db::start")
        build_bowtie2_db(bt2db_dir, bt2db_name, contigs_files, args.num_cores)
        tsprint(f"MIDAS::build_bowtie2db::finish")

        tsprint(f"MIDAS::bowtie2_align::start")
        repgenome_bamfile = sample.get_target_layout("snps_repgenomes_bam")
        bowtie2_align(bt2db_dir, bt2db_name, repgenome_bamfile, args)
        samtools_index(repgenome_bamfile, args.debug, args.num_cores)
        tsprint(f"MIDAS::bowtie2_align::finish")

        tsprint(f"MIDAS::filter_bam::start")
        args_list = []
        for species_id in species_ids_of_interest:
            repgenome_bamfile = sample.get_target_layout("snps_repgenomes_bam") # input
            filtered_bamfile = sample.get_target_layout("species_bam", species_id) # <species_id> output
            args_list.append((species_id, repgenome_bamfile, filtered_bamfile))
        list_of_contig_aln_stats = multiprocessing_map(filter_bam, args_list, args.num_cores)

        cores_per_species = max(4, floor(args.num_cores/species_counts))
        multiprocessing_map(sort_bam, [(species_id, cores_per_species) for species_id in species_ids_of_interest], floor(args.num_cores/cores_per_species))
        tsprint(f"MIDAS::filter_bam::finish")

        tsprint(f"MIDAS::multiprocessing_map::start")
        chunks_pileup_summary = multiprocessing_map(process_chunk_of_sites, arguments_list, args.num_cores)
        tsprint(f"MIDAS::multiprocessing_map::finish")

        tsprint(f"MIDAS::write_species_pileup_summary::start")
        snps_summary_fp = sample.get_target_layout("snps_summary")
        snps_chunk_summary_fp = sample.get_target_layout("snps_chunk_summary")

        dict_of_chunk_aln_stats = compute_chunk_aln_summary(list_of_contig_aln_stats, species_ids_of_interest)
        write_species_pileup_summary(chunks_pileup_summary, snps_summary_fp, snps_chunk_summary_fp, dict_of_chunk_aln_stats)
        tsprint(f"MIDAS::write_species_pileup_summary::finish")

    except Exception as error:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error. Specify --debug flag to keep.")
            sample.remove_dirs(["outdir", "tempdir"])
            if not args.prebuilt_bowtie2_indexes:
                sample.remove_dirs(["bt2_indexes_dir"])
        raise error


@register_args
def main(args):
    tsprint(f"Single nucleotide polymorphisms calling in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    run_snps(args)
