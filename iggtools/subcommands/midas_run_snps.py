#!/usr/bin/env python3
import json
import os
import multiprocessing

from collections import defaultdict
from pysam import AlignmentFile  # pylint: disable=no-name-in-module

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, num_physical_cores, InputStream, OutputStream, multiprocessing_map, command, cat_files, select_from_tsv, multithreading_map
from iggtools.models.uhgg import MIDAS_IGGDB
from iggtools.common.bowtie2 import build_bowtie2_db, bowtie2_align, samtools_index, bowtie2_index_exists, _keep_read
from iggtools.params.schemas import snps_profile_schema, snps_pileup_schema, format_data, snps_chunk_summary_schema
from iggtools.models.sample import Sample
from iggtools.models.species import Species


DEFAULT_MARKER_DEPTH = 5.0
DEFAULT_ALN_MAPID = 94.0
DEFAULT_ALN_MAPQ = 20
DEFAULT_ALN_READQ = 20
DEFAULT_ALN_BASEQ = 30
DEFAULT_ALN_COV = 0.75
DEFAULT_ALN_TRIM = 0
DEFAULT_CHUNK_SIZE = 50000


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
    subparser.add_argument('--midas_iggdb',
                           dest='midas_iggdb',
                           type=str,
                           metavar="CHAR",
                           help=f"local MIDAS DB which mirrors the s3 IGG db")

    # Species related
    subparser.add_argument('--marker_depth',
                           type=float,
                           dest='marker_depth',
                           metavar='FLOAT',
                           default=DEFAULT_MARKER_DEPTH,
                           help=f"Include species with >X coverage ({DEFAULT_MARKER_DEPTH})")
    subparser.add_argument('--species_list',
                           dest='species_list',
                           type=str,
                           metavar="CHAR",
                           help=f"Comma separated list of species ids")

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
    # File related
    subparser.add_argument('--sparse',
                           action='store_true',
                           default=True,
                           help=f"Omit zero rows from output.")
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
    return main_func


def keep_read(aln):
    global global_args
    args = global_args
    return _keep_read(aln, args.aln_mapid, args.aln_readq, args.aln_mapq, args.aln_cov)


def design_chunks_per_species(args):
    sp, midas_iggdb, chunk_size = args
    return sp.design_snps_chunks(midas_iggdb, chunk_size)


def design_chunks(species_ids_of_interest, midas_iggdb, chunk_size):
    """ Chunks of continuous genomics sites, indexed by species_id, chunk_id """

    global semaphore_for_species
    global dict_of_species

    # Read-only global variables
    semaphore_for_species = dict()
    dict_of_species = {species_id: Species(species_id) for species_id in species_ids_of_interest}

    # Design chunks structure per species
    flags = multithreading_map(design_chunks_per_species, [(sp, midas_iggdb, chunk_size) for sp in dict_of_species.values()], 4)
    #flags = [sp.design_snps_chunks(midas_iggdb, chunk_size) for sp in dict_of_species.values()]
    assert all(flags)

    # Prioritize chunks need to compute read counts across species
    # TODO: if we don't need to priority across species, then can put the following part into design_snps_chunks
    pri_list = []
    reg_list = []
    for species_id, sp in dict_of_species.items():
        priority_chunks = sp.priority_chunks
        num_of_sites_chunks = sp.num_of_sites_chunks
        for chunk_id in range(0, num_of_sites_chunks):
            if chunk_id in priority_chunks:
                pri_list.append((species_id, chunk_id))
            else:
                reg_list.append((species_id, chunk_id))
        reg_list.append((species_id, -1))

        # Create a semaphore with number_of_chunks of elements
        semaphore_for_species[species_id] = multiprocessing.Semaphore(num_of_sites_chunks)
        for _ in range(num_of_sites_chunks):
            semaphore_for_species[species_id].acquire()

    arguments_list = pri_list + reg_list
    return arguments_list


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

    try:
        species_id, chunk_id = packed_args
        sp = dict_of_species[species_id]
        ret = []
        for pargs in sp.chunks_of_sites[chunk_id]:
            ret.append(pileup_per_unit(pargs))
        return ret
    finally:
        semaphore_for_species[species_id].release() # no deadlock


def pileup_per_unit(packed_args):
    """ Pileup for continuous of one contig in one chunk """

    global global_args
    global dict_of_species
    global sample

    # [contig_start, contig_end)
    species_id, chunk_id, contig_id, contig_start, contig_end, count_flag = packed_args #contig_seq

    repgenome_bamfile = sample.get_target_layout("snps_repgenomes_bam")
    headerless_sliced_path = sample.get_target_layout("chunk_pileup", species_id, chunk_id)
    contig_seq = dict_of_species[species_id].contigs[contig_id]["seq"]

    zero_rows_allowed = not global_args.sparse
    current_chunk_size = contig_end - contig_start
    aligned_reads = 0
    mapped_reads = 0

    with AlignmentFile(repgenome_bamfile) as bamfile:
        counts = bamfile.count_coverage(contig_id, contig_start, contig_end,
                                        quality_threshold=global_args.aln_baseq, # min_quality_threshold a base has to reach to be counted.
                                        read_callback=keep_read) # select a call-back to ignore reads when counting
        if count_flag:
            # Single read could cover the chunk boundaries, and to avoid overcounting of boundary reads,
            # we only compute the aligned_reads per contig once
            aligned_reads = bamfile.count(contig_id)
            mapped_reads = bamfile.count(contig_id, read_callback=keep_read)

    # aln_stats need to be passed from child process back to parents
    aln_stats = {
        "species_id": species_id,
        "chunk_id": chunk_id,
        "contig_id": contig_id,
        "chunk_length": current_chunk_size,
        "aligned_reads": aligned_reads,
        "mapped_reads": mapped_reads,
        "contig_total_depth": 0,
        "contig_covered_bases": 0
    }

    # TODO: instead of write to file, save into memory.
    with open(headerless_sliced_path, "a") as stream:
        for within_chunk_index in range(0, current_chunk_size):
            depth = sum([counts[nt][within_chunk_index] for nt in range(4)])
            count_a = counts[0][within_chunk_index]
            count_c = counts[1][within_chunk_index]
            count_g = counts[2][within_chunk_index]
            count_t = counts[3][within_chunk_index]

            ref_pos = within_chunk_index + contig_start
            ref_allele = contig_seq[ref_pos]
            row = (contig_id, ref_pos + 1, ref_allele, depth, count_a, count_c, count_g, count_t)

            aln_stats["contig_total_depth"] += depth
            if depth > 0:
                aln_stats["contig_covered_bases"] += 1
            if depth > 0 or zero_rows_allowed:
                stream.write("\t".join(map(format_data, row)) + "\n")
        assert within_chunk_index+contig_start == contig_end-1, f"compute_pileup_per_chunk::index mismatch error for {contig_id}."
    return aln_stats


def merge_chunks_per_species(species_id):
    """ merge the pileup results from chunks into one file per species """

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
    if False: #not global_args.debug:
        tsprint(f"Deleting temporary sliced pileup files for {species_id}.")
        for s_file in list_of_chunks_pileup:
            command(f"rm -rf {s_file}", quiet=True)

    # return a status flag
    # the path should be computable somewhere else
    return True


def write_species_pileup_summary(chunks_pileup_summary, outfile, chunk_output):
    """ Collect species pileup aln stats from all chunks and write to file """

    species_pileup_summary = defaultdict(dict)
    with OutputStream(chunk_output) as stream:
        stream.write("\t".join(snps_chunk_summary_schema.keys()) + "\n")
        for records in chunks_pileup_summary:
            if records is True:
                continue
            for record in records:
                stream.write("\t".join(map(format_data, record.values())) + "\n")
                species_id = record["species_id"]
                if species_id not in species_pileup_summary:
                    species_pileup_summary[species_id] = {
                        "species_id": species_id,
                        "genome_length": 0,
                        "covered_bases": 0,
                        "total_depth": 0,
                        "aligned_reads":0,
                        "mapped_reads":0,
                        "fraction_covered": 0.0,
                        "mean_coverage": 0.0
                        }

                curr_species_pileup = species_pileup_summary.get(species_id)
                curr_species_pileup["genome_length"] += record["chunk_length"]
                curr_species_pileup["total_depth"] += record["contig_total_depth"]
                curr_species_pileup["covered_bases"] += record["contig_covered_bases"]
                curr_species_pileup["aligned_reads"] += record["aligned_reads"]
                curr_species_pileup["mapped_reads"] += record["mapped_reads"]

    # Secondary round compute: need to loop over all chunks of sites to compute fraction_covered
    for species_id in species_pileup_summary.keys():
        curr_species_pileup = species_pileup_summary.get(species_id)
        if curr_species_pileup["genome_length"] > 0:
            curr_species_pileup["fraction_covered"] = curr_species_pileup["covered_bases"] / curr_species_pileup["genome_length"]
        if curr_species_pileup["covered_bases"] > 0:
            curr_species_pileup["mean_coverage"] = curr_species_pileup["total_depth"] / curr_species_pileup["covered_bases"]

    # Write to file
    with OutputStream(outfile) as stream:
        stream.write("\t".join(snps_profile_schema.keys()) + "\n")
        for record in species_pileup_summary.values():
            stream.write("\t".join(map(format_data, record.values())) + "\n")


def midas_run_snps(args):

    try:
        global sample
        sample = Sample(args.sample_name, args.midas_outdir, "snps")
        sample.create_dirs(["outdir", "tempdir"], args.debug, quiet=True)


        global global_args
        global_args = args


        species_list = args.species_list.split(",") if args.species_list else []
        if args.prebuilt_bowtie2_indexes:
            bt2_db_dir = os.path.dirname(args.prebuilt_bowtie2_indexes)
            bt2_db_name = os.path.basename(args.prebuilt_bowtie2_indexes)
            assert bowtie2_index_exists(bt2_db_dir, bt2_db_name), f"Provided {bt2_db_dir}/{bt2_db_name} don't exist."

            # We only need a list of species that we need to pull the
            assert (args.prebuilt_bowtie2_species and os.path.exists(args.prebuilt_bowtie2_species)), f"Need to provide list of speices used to build the provided Bowtie2 indexes."
            tsprint(f"Read in list of species used to build provided bowtie2 indexes {bt2_db_dir}/{bt2_db_name}")
            bt2_species_list = []
            with InputStream(args.prebuilt_bowtie2_species) as stream:
                for species_id in select_from_tsv(stream, schema={"species_id": str}):
                    bt2_species_list.append(species_id[0])

            # Update species_list: either particular species of interest or species in the bowtie2 indexes
            species_list = list(set(species_list) & set(bt2_species_list)) if species_list else bt2_species_list
            # Note: the index was only for the purpose of same bowtie2 index. but the species_ids_of_interest per sample
            # can still be based on sample itself.
            # We also don't want too many empty species in our parsing stage.
        else:
            sample.create_dirs(["bt2_indexes_dir"], args.debug, quiet=True)
            bt2_db_dir = sample.get_target_layout("bt2_indexes_dir")
            bt2_db_name = "repgenomes"


        # Select abundant species present in the sample for SNPs calling
        species_ids_of_interest = species_list if args.marker_depth == -1 else sample.select_species(args.marker_depth, species_list)
        species_counts = len(species_ids_of_interest)
        sample.create_species_subdirs(species_ids_of_interest, "temp", args.debug, quiet=True)
        assert species_counts > 0, f"No (specified) species pass the marker_depth filter, please adjust the marker_depth or species_list"
        tsprint(species_ids_of_interest)


        # Fetch representative genome fastas for each species (multiprocessing)
        tsprint(f"CZ::design_chunks::start")
        num_cores = min(args.num_cores, species_counts)
        midas_iggdb = MIDAS_IGGDB(args.midas_iggdb if args.midas_iggdb else sample.get_target_layout("midas_iggdb_dir"), num_cores)
        arguments_list = design_chunks(species_ids_of_interest, midas_iggdb, args.chunk_size)
        tsprint(f"CZ::design_chunks::finish")


        # Build Bowtie indexes for species in the restricted species profile
        contigs_files = midas_iggdb.fetch_files("prokka_genome", species_ids_of_interest)
        tsprint(f"CZ::build_bowtie2_indexes::start")
        build_bowtie2_db(bt2_db_dir, bt2_db_name, contigs_files, args.num_cores)
        tsprint(f"CZ::build_bowtie2_indexes::finish")
        # Perhaps avoid this giant conglomerated file, fetching instead submaps for each species.
        # TODO: Also colocate/cache/download in master for multiple slave subcommand invocations


        tsprint(f"CZ::bowtie2_align::start")
        repgenome_bamfile = sample.get_target_layout("snps_repgenomes_bam")
        bowtie2_align(bt2_db_dir, bt2_db_name, repgenome_bamfile, args)
        samtools_index(repgenome_bamfile, args.debug, args.num_cores)
        tsprint(f"CZ::bowtie2_align::finish")


        tsprint(f"CZ::multiprocessing_map::start")
        chunks_pileup_summary = multiprocessing_map(process_one_chunk_of_sites, arguments_list, args.num_cores)
        tsprint(f"CZ::multiprocessing_map::finish")


        tsprint(f"CZ::write_species_pileup_summary::start")
        write_species_pileup_summary(chunks_pileup_summary, sample.get_target_layout("snps_summary"), sample.get_target_layout("snps_chunk_summary"))
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
