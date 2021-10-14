#!/usr/bin/env python3
import json
import os
import multiprocessing
from operator import itemgetter
from collections import defaultdict

import numpy as np
from pysam import AlignmentFile  # pylint: disable=no-name-in-module

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, num_physical_cores, InputStream, OutputStream, multiprocessing_map, command, cat_files, select_from_tsv, multithreading_map
from iggtools.models.midasdb import MIDAS_DB
from iggtools.common.bowtie2 import build_bowtie2_db, bowtie2_align, samtools_sort, samtools_index, bowtie2_index_exists, _keep_read
from iggtools.params.schemas import snps_profile_schema, snps_pileup_schema, format_data, snps_chunk_summary_schema
from iggtools.models.sample import Sample
from iggtools.models.species import Species, collect_units_per_chunk, parse_species
from iggtools.common.snvs import call_alleles, reference_overlap, update_overlap, mismatches_within_overlaps


DEFAULT_MARKER_DEPTH = 5.0
DEFAULT_MARKER_MEDIAN_DEPTH = 0

DEFAULT_ALN_MAPID = 94.0
DEFAULT_ALN_MAPQ = 10
DEFAULT_ALN_COV = 0.75
DEFAULT_ALN_READQ = 20
DEFAULT_ALN_BASEQ = 30
DEFAULT_ALN_TRIM = 0

DEFAULT_CHUNK_SIZE = 50000
DEFAULT_MAX_FRAGLEN = 5000

DEFAULT_SITE_DEPTH = 2
DEFAULT_SNP_MAF = 0.05


def register_args(main_func):
    subparser = add_subcommand('midas_run_snps', main_func, help='Predict single-nucleotide-polymorphism')

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
    subparser.add_argument('--midas_db',
                           dest='midas_db',
                           type=str,
                           metavar="CHAR",
                           help=f"local MIDAS DB which mirrors the s3 IGG db")

    # Species related
    subparser.add_argument('--species_list',
                           dest='species_list',
                           type=str,
                           metavar="CHAR",
                           help=f"Comma separated list of species ids")
    subparser.add_argument('--select_by',
                           dest='select_by',
                           type=str,
                           metavar="CHAR",
                           default="median_marker_coverage",
                           choices=['median_marker_coverage', 'marker_coverage', 'unique_fraction_covered', "marker_relative_abundance"],
                           help=f"Column from species_profile based on which to select species.")
    subparser.add_argument('--select_threshold',
                           dest='select_threshold',
                           type=float,
                           metavar="FLOAT",
                           default=DEFAULT_MARKER_MEDIAN_DEPTH,
                           help=f"Include species with > X median SGC (median) marker coverage ({DEFAULT_MARKER_MEDIAN_DEPTH})")

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
                           help=f"Maximum fragment length for paired reads.")


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
                           default=num_physical_cores,
                           help=f"Number of physical cores to use ({num_physical_cores})")

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


def design_chunks_per_species(args):
    sp, midas_db, chunk_size = args
    return sp.design_snps_chunks(midas_db, chunk_size)


def design_chunks(species_ids_of_interest, midas_db, chunk_size):
    """ Chunks of continuous genomics sites, indexed by species_id, chunk_id """

    global semaphore_for_species
    global dict_of_species

    # Read-only global variables
    semaphore_for_species = dict()
    dict_of_species = {species_id: Species(species_id) for species_id in species_ids_of_interest}

    # Design chunks structure per species
    flags = multithreading_map(design_chunks_per_species, [(sp, midas_db, chunk_size) for sp in dict_of_species.values()], 4) #<---
    assert all(flags)

    # Sort species by the largest contig length
    sorted_tuples_of_species = sorted(((sp.id, sp.max_contig_length) for sp in dict_of_species.values()), key=itemgetter(1), reverse=True)

    arguments_list = []
    for species_id, _ in sorted_tuples_of_species:
        sp = dict_of_species[species_id]

        # The structure of the chunks depends on the representative genome sequences
        num_of_sites_chunks = sp.num_of_sites_chunks
        for chunk_id in range(0, num_of_sites_chunks):
            arguments_list.append((species_id, chunk_id))

        # Create a semaphore with number_of_chunks of elements
        semaphore_for_species[species_id] = multiprocessing.Semaphore(num_of_sites_chunks)
        for _ in range(num_of_sites_chunks):
            semaphore_for_species[species_id].acquire()

    for species_id in dict_of_species.keys():
        arguments_list.append((species_id, -1))

    return arguments_list


def filter_bam_by_single_read(pargs):
    """ Filter given BAM file with propely paired reads for given species and write to file """

    global  global_args
    global dict_of_species
    global sample

    repbamfile, species_id = pargs

    tsprint(f"  CZ::filter_bam_by_single_read::{species_id}-0::start filter_bam_by_single_read")

    # List of contigs for given species
    list_of_contig_ids = list(dict_of_species[species_id].contigs.keys())

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
            mapped_reads = 0
            for aln in infile.fetch(contig_id):
                aligned_reads += 1
                if keep_read(aln):
                    mapped_reads += 1
                    read = "1" if aln.is_read1 else "2"
                    filtered_alns_dict[f"{aln.query_name}_{read}"] = aln

            reads_stats["aligned_reads"][contig_id] = aligned_reads
            reads_stats["mapped_reads"][contig_id] = mapped_reads

    # Write filtered alignments to file
    template_bam = AlignmentFile(repbamfile, "rb")
    filtered_bam = AlignmentFile(sample.get_target_layout("species_bam", species_id), "wb", template=template_bam)
    for aln in filtered_alns_dict.values():
        filtered_bam.write(aln)
    filtered_bam.close()

    tsprint(f"  CZ::filter_bam_by_single_read::{species_id}-0::finish filter_bam_by_single_read")
    return reads_stats


def filter_bam_by_proper_pair(pargs):
    """ Filter given BAM file with propely paired reads for given species and write to file """

    global global_args
    global dict_of_species
    global sample

    repbamfile, species_id = pargs

    tsprint(f"  CZ::filter_bam_by_proper_pair::{species_id}-0::start filter_bam_by_proper_pair")

    # List of contigs for given species
    list_of_contig_ids = list(dict_of_species[species_id].contigs.keys())

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

                overlap_pass = True
                if reads_overlap:
                    # Keep the FWD read, split the REV reads
                    (nm_out_rev, nm_in_rev, ngaps_ri_rev, ngaps_ro_rev) = mismatches_within_overlaps(alns["rev"], reads_overlap, "rev")
                    #assert nm_out_rev + nm_in_rev == dict(alns["rev"].tags)['NM']

                    # Keep the REV read, split the FWD reads
                    (nm_out_fwd, nm_in_fwd, ngaps_ri_fwd, ngaps_ro_fwd) = mismatches_within_overlaps(alns["fwd"], reads_overlap, "fwd")
                    #assert nm_out_fwd + nm_in_fwd == dict(alns["fwd"].tags)['NM']

                    # Update the overlap by substracting the number of gaps in the fwd overlap region
                    reads_overlap = reads_overlap - ngaps_ri_fwd

                    # For repeats regions, paired-end reads can be aligned with many gaps, resulting in high mismatches within the overlapping region
                    # Only keep aligned pairs indicating from the same DNA fragment
                    if abs(nm_in_fwd - nm_in_rev) > 1:
                        overlap_pass = False
                        continue #<-----------

                    mismatches = dict(alns["fwd"].tags)['NM'] + nm_out_rev

                    # Update the aligned_length to compute the mapid
                    align_len = alns["rev"].query_alignment_length + alns["fwd"].query_alignment_length - reads_overlap
                    align_len_no_gaps = align_len - ngaps_ro_rev - ngaps_ro_fwd - ngaps_ri_fwd

                    # To avoid overcounting site depth for the overlapping region,
                    # "The higher quality base is used and the lower-quality base is set to BQ=0."
                    b1 = alns["fwd"].query_alignment_end - reads_overlap
                    b2 = alns["rev"].query_alignment_start + reads_overlap - 1

                    # TODO: loose end => this gives me error using ibd data
                    debug_string = "\t".join([str(b1), str(b2), str(reads_overlap), str(len(alns["fwd"].query_alignment_sequence[b1:])), str(len(alns["rev"].query_alignment_sequence[:b2+1])), str(overlap_pass)])
                    #assert reads_overlap == len(alns["fwd"].query_alignment_sequence[b1:]), debug_string
                    #assert reads_overlap == len(alns["rev"].query_alignment_sequence[:b2+1]), debug_string
                    #assert len(alns["fwd"].query_alignment_sequence[b1:]) == len(alns["rev"].query_alignment_sequence[:b2+1])
                    # also here
                    #if False and reads_overlap != len(alns["rev"].query_alignment_sequence[:b2+1]) and overlap_pass:
                    #    debug_overlap(alns)
                    #    print(debug_string)

                    # Only use the higher quality base in the overlap region for downstream pileup
                    f = alns["fwd"].query_qualities[b1:]
                    r = alns["rev"].query_qualities[:b2+1]
                    for i, _ in enumerate(zip(f, r)):
                        (x, y) = _
                        if x >= y:
                            r[i] = 0
                        else:
                            f[i] = 0
                    alns["fwd"].query_qualities[b1:] = f
                    alns["rev"].query_qualities[:b2+1] = r

                    mapid = 100 * (align_len - mismatches) / float(align_len)
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
    filtered_bam = AlignmentFile(sample.get_target_layout("species_bam", species_id), "wb", template=template_bam)
    for query_name, alns in filtered_alns_dict.items():
        filtered_bam.write(alns["fwd"])
        filtered_bam.write(alns["rev"])
    filtered_bam.close()

    tsprint(f"  CZ::filter_bam_by_proper_pair::{species_id}-0::finish filter_bam_by_proper_pair")
    return reads_stats


def process_one_chunk_of_sites(packed_args):
    """ Process one chunk: either pileup or merge and write results to disk """

    species_id, chunk_id = packed_args

    if chunk_id == -1:
        global semaphore_for_species
        global dict_of_species
        sp = dict_of_species[species_id]

        tsprint(f"  CZ::process_one_chunk_of_sites::{species_id}-{chunk_id}::wait merge_chunks_per_species")
        for _ in range(sp.num_of_sites_chunks):
            semaphore_for_species[species_id].acquire()
        tsprint(f"  CZ::process_one_chunk_of_sites::{species_id}-{chunk_id}::start merge_chunks_per_species")
        ret = merge_chunks_per_species(species_id)
        tsprint(f"  CZ::process_one_chunk_of_sites::{species_id}-{chunk_id}::finish merge_chunks_per_species")
        return ret

    tsprint(f"  CZ::process_one_chunk_of_sites::{species_id}-{chunk_id}::start compute_pileup_per_chunk")
    ret = compute_pileup_per_chunk(packed_args)
    tsprint(f"  CZ::process_one_chunk_of_sites::{species_id}-{chunk_id}::finish compute_pileup_per_chunk")

    return ret


def compute_pileup_per_chunk(packed_args):
    """ Pileup for one chunk, potentially contain multiple contigs """

    global semaphore_for_species
    global dict_of_species
    global sample
    global global_args

    try:
        species_id, chunk_id = packed_args
        sp = dict_of_species[species_id]
        ret = []
        for pargs in sp.chunks_of_sites[chunk_id]:
            ret.append(pileup_per_unit(pargs))

        contig_counts_per_chunk = len(sp.chunks_of_sites[chunk_id])
        collect_units_per_chunk(sample, contig_counts_per_chunk, species_id, chunk_id, "chunk_pileup")

        return ret
    finally:
        semaphore_for_species[species_id].release() # no deadlock


def pileup_per_unit(packed_args):
    """ Pileup for continuous positions of one contig in one chunk """

    global global_args
    global dict_of_species
    global sample

    # [contig_start, contig_end)
    species_id, chunk_id, contig_id, contig_start, contig_end, count_flag, within_chunk_cid = packed_args

    repgenome_bamfile = sample.get_target_layout("species_sorted_bam", species_id)
    headerless_sliced_path = sample.get_target_layout("chunk_pileup_perc", species_id, chunk_id, within_chunk_cid)

    current_chunk_size = contig_end - contig_start
    contig_seq = dict_of_species[species_id].contigs[contig_id]["seq"]

    with AlignmentFile(repgenome_bamfile) as bamfile:
        # min_quality_threshold a base has to reach to be counted.
        counts = bamfile.count_coverage(contig_id, contig_start, contig_end, quality_threshold=global_args.aln_baseq)

    # aln_stats need to be passed from child process back to parents
    aln_stats = {
        "species_id": species_id,
        "chunk_id": chunk_id,
        "contig_id": contig_id,
        "chunk_length": current_chunk_size,
        "aligned_reads": 0, #aligned_reads if count_flag and not global_args.paired_only else 0,
        "mapped_reads": 0, #mapped_reads if count_flag and not global_args.paired_only else 0,
        "contig_total_depth": 0,
        "contig_covered_bases": 0,
    }

    with OutputStream(headerless_sliced_path) as stream:
        for within_chunk_index in range(0, current_chunk_size):
            depth = sum([counts[nt][within_chunk_index] for nt in range(4)])
            count_a = counts[0][within_chunk_index]
            count_c = counts[1][within_chunk_index]
            count_g = counts[2][within_chunk_index]
            count_t = counts[3][within_chunk_index]

            ref_pos = within_chunk_index + contig_start
            ref_allele = contig_seq[ref_pos]
            row = [contig_id, ref_pos + 1, ref_allele, depth, count_a, count_c, count_g, count_t]

            aln_stats["contig_total_depth"] += depth
            if depth < global_args.site_depth:
                continue
            aln_stats["contig_covered_bases"] += 1

            rc_ACGT = [count_a, count_c, count_g, count_t]
            tuple_of_alleles = zip(['A', 'C', 'G', 'T'], rc_ACGT)
            major_allele, minor_allele, snp_type, allele_counts = call_alleles(tuple_of_alleles, depth, global_args.snp_maf)

            if allele_counts == 0:
                continue
            
            major_index = 'ACGT'.index(major_allele)
            minor_index = 'ACGT'.index(minor_allele)
            major_allelefreq = rc_ACGT[major_index] / depth
            minor_allelefreq = 0.0 if major_index == minor_index else rc_ACGT[minor_index] / depth

            row.extend([major_allele, minor_allele, major_allelefreq, minor_allelefreq, allele_counts])

            stream.write("\t".join(map(format_data, row)) + "\n")
        assert within_chunk_index+contig_start == contig_end-1, f"compute_pileup_per_chunk::index mismatch error for {contig_id}."

    # Delete temporary bam file <== TODO
    if global_args.paired_only and False:
        command(f"rm -rf {repgenome_bamfile}", quiet=False)
        command(f"rm -rf {repgenome_bamfile}.bai", quiet=True)

    return aln_stats


def merge_chunks_per_species(species_id):
    """ merge the pileup results from chunks into one file per species """

    global global_args
    global sample
    global dict_of_species

    sp = dict_of_species[species_id]
    number_of_chunks = sp.num_of_sites_chunks

    list_of_chunks_pileup = [sample.get_target_layout("chunk_pileup", species_id, chunk_id) for chunk_id in range(0, number_of_chunks)]
    species_snps_pileup_file = sample.get_target_layout("snps_pileup", species_id)

    with OutputStream(species_snps_pileup_file) as stream:
        stream.write('\t'.join(snps_pileup_schema.keys()) + '\n')
    cat_files(list_of_chunks_pileup, species_snps_pileup_file, 20)

    # The chunk_pilup_path will be used in merge_midas_snps.
    if not global_args.debug:
        tsprint(f"Deleting temporary sliced pileup files for {species_id}.")
        for s_file in list_of_chunks_pileup:
            command(f"rm -rf {s_file}", quiet=True)

    # return a status flag
    # the path should be computable somewhere else
    return True


def assign_contig_reads_to_chunks(lalns_stats_by_contig, species_ids_of_interest):
    global dict_of_species
    dchunk_alns_stats = dict()

    for spidx, species_id in enumerate(species_ids_of_interest):
        sp = dict_of_species[species_id]

        cc_to_ch = defaultdict(lambda: defaultdict(dict))
        for chunk_id, tchunks_list in sp.chunks_of_sites.items():
            if chunk_id == -1:
                continue
            #[(species_id, chunk_id, contig_id, ci, ci+chunk_size, count_flag, 0)]
            for loc in tchunks_list:
                contig_id = loc[2]
                count_flag = loc[5]

                aligned_reads = lalns_stats_by_contig[spidx]["aligned_reads"][contig_id]
                mapped_reads = lalns_stats_by_contig[spidx]["mapped_reads"][contig_id]

                if contig_id not in cc_to_ch[chunk_id]:
                    cc_to_ch[chunk_id][contig_id] = {"aligned_reads": 0, "mapped_reads": 0}

                cc_to_ch[chunk_id][contig_id]["aligned_reads"] = aligned_reads if count_flag else 0
                cc_to_ch[chunk_id][contig_id]["mapped_reads"] = mapped_reads if count_flag else 0

        dchunk_alns_stats[species_id] = cc_to_ch
    return dchunk_alns_stats


def write_species_pileup_summary(chunks_pileup_summary, snps_summary_outfile, chunk_output, dchunk_alns_stats):
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

                record["aligned_reads"] = dchunk_alns_stats[species_id][chunk_id][contig_id]["aligned_reads"]
                record["mapped_reads"] = dchunk_alns_stats[species_id][chunk_id][contig_id]["mapped_reads"]

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


def midas_run_snps(args):

    try:
        global global_args
        global_args = args

        global sample
        sample = Sample(args.sample_name, args.midas_outdir, "snps")
        sample.create_dirs(["outdir", "tempdir"], args.debug, quiet=True)

        species_list = parse_species(args)

        # Prepare Bowtie2 genome database path and name, prebuilt or to be built.
        if args.prebuilt_bowtie2_indexes:
            bt2db_dir = os.path.dirname(args.prebuilt_bowtie2_indexes)
            bt2db_name = os.path.basename(args.prebuilt_bowtie2_indexes)
            assert bowtie2_index_exists(bt2db_dir, bt2db_name), f"Provided {bt2db_dir}/{bt2db_name} don't exist."

            # Required list of species used to build bowtie2 database to fetch the genome sequences.
            assert (args.prebuilt_bowtie2_species and os.path.exists(args.prebuilt_bowtie2_species)), f"Need to provide list of speices used to build the provided Bowtie2 indexes."
            tsprint(f"Read in list of species used to build provided bowtie2 indexes {bt2db_dir}/{bt2db_name}")
            bt2_species_list = []
            with InputStream(args.prebuilt_bowtie2_species) as stream:
                for species_id in select_from_tsv(stream, schema={"species_id": str}):
                    bt2_species_list.append(species_id[0])

            # Update the species list
            species_list = list(set(species_list) & set(bt2_species_list)) if species_list else bt2_species_list
        else:
            sample.create_dirs(["bt2_indexes_dir"], args.debug, quiet=True)
            bt2db_dir = sample.get_target_layout("bt2_indexes_dir")
            bt2db_name = "repgenomes"


        # Pileup species: abundant and/or listed species. We don't recommend pileup on too many empty species.
        species_ids_of_interest = species_list if args.select_threshold == -1 else sample.select_species(args, species_list)

        species_counts = len(species_ids_of_interest)

        sample.create_species_subdirs(species_ids_of_interest, "temp", args.debug, quiet=True)
        assert species_counts > 0, f"No (specified) species pass the marker_depth filter, please adjust the marker_depth or species_list"
        tsprint(len(species_ids_of_interest))


        # Fetch representative genome fastas for each species (multiprocessing)
        tsprint(f"CZ::design_chunks::start")
        num_cores = min(args.num_cores, species_counts)
        midas_db = MIDAS_DB(args.midas_db if args.midas_db else sample.get_target_layout("midas_db_dir")) #<---- all species
        arguments_list = design_chunks(species_ids_of_interest, midas_db, args.chunk_size)
        tsprint(f"CZ::design_chunks::finish")


        # Build Bowtie indexes for species in the restricted species profile
        contigs_files = midas_db.fetch_files("prokka_genome", species_ids_of_interest)
        tsprint(f"CZ::build_bowtie2_indexes::start")
        build_bowtie2_db(bt2db_dir, bt2db_name, contigs_files, args.num_cores)
        tsprint(f"CZ::build_bowtie2_indexes::finish")


        tsprint(f"CZ::bowtie2_align::start")
        repgenome_bamfile = sample.get_target_layout("snps_repgenomes_bam")
        bowtie2_align(bt2db_dir, bt2db_name, repgenome_bamfile, args)
        samtools_index(repgenome_bamfile, args.debug, args.num_cores)
        tsprint(f"CZ::bowtie2_align::finish")


        tsprint(f"CZ::filter_bam_by_proper_pair::start")
        args_list = [(repgenome_bamfile, species_id) for species_id in species_ids_of_interest]
        if args.paired_only:
            lalns_stats_by_contig = multiprocessing_map(filter_bam_by_proper_pair, args_list, args.num_cores)
        else:
            lalns_stats_by_contig = multiprocessing_map(filter_bam_by_single_read, args_list, args.num_cores)

        for species_id in species_ids_of_interest:
            species_bam = sample.get_target_layout("species_bam", species_id)
            species_sorted_bam = sample.get_target_layout("species_sorted_bam", species_id)
            samtools_sort(species_bam, species_sorted_bam, args.debug, args.num_cores)
            samtools_index(species_sorted_bam, args.debug, args.num_cores)
            command(f"rm -rf {species_bam}", quiet=True)
        tsprint(f"CZ::filter_bam_by_proper_pair::finish")


        tsprint(f"CZ::multiprocessing_map::start")
        chunks_pileup_summary = multiprocessing_map(process_one_chunk_of_sites, arguments_list, args.num_cores)
        tsprint(f"CZ::multiprocessing_map::finish")


        tsprint(f"CZ::write_species_pileup_summary::start")
        snps_summary_fp = sample.get_target_layout("snps_summary")
        snps_chunk_summary_fp = sample.get_target_layout("snps_chunk_summary")
        dchunk_alns_stats = assign_contig_reads_to_chunks(lalns_stats_by_contig, species_ids_of_interest)
        write_species_pileup_summary(chunks_pileup_summary, snps_summary_fp, snps_chunk_summary_fp, dchunk_alns_stats)
        tsprint(f"CZ::write_species_pileup_summary::finish")


    except Exception as error:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error. Specify --debug flag to keep.")
            sample.remove_dirs(["outdir", "tempdir"])
            if not args.prebuilt_bowtie2_indexes:
                sample.remove_dirs(["bt2_indexes_dir"])
        raise error


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_run_snps(args)
