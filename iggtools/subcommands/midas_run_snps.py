import json
import os
import multiprocessing

from collections import defaultdict
from math import ceil
import numpy as np
from pysam import AlignmentFile  # pylint: disable=no-name-in-module
import Bio.SeqIO

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, num_physical_cores, InputStream, OutputStream, multiprocessing_map, download_reference, command, cat_files, select_from_tsv
from iggtools.params import outputs
from iggtools.models.uhgg import UHGG
from iggtools.common.bowtie2 import build_bowtie2_db, bowtie2_align, samtools_index, bowtie2_index_exists
from iggtools.params.schemas import snps_profile_schema, snps_pileup_schema, format_data
from iggtools.models.sample import Sample


DEFAULT_ALN_COV = 0.75
DEFAULT_GENOME_COVERAGE = 3.0
DEFAULT_ALN_MAPID = 94.0
DEFAULT_ALN_MAPQ = 20
DEFAULT_ALN_READQ = 20
DEFAULT_ALN_BASEQ = 30
DEFAULT_ALN_TRIM = 0
DEFAULT_CHUNK_SIZE = 40000


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
    # Species related
    subparser.add_argument('--genome_coverage',
                           type=float,
                           dest='genome_coverage',
                           metavar='FLOAT',
                           default=DEFAULT_GENOME_COVERAGE,
                           help=f"Include species with >X coverage ({DEFAULT_GENOME_COVERAGE})")
    subparser.add_argument('--species_list',
                           dest='species_list',
                           type=str,
                           metavar="CHAR",
                           help=f"Comma separated list of species ids")
    subparser.add_argument('--local_bowtie2_indexes',
                           dest='local_bowtie2_indexes',
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
    # File related
    subparser.add_argument('--sparse',
                           action='store_true',
                           default=False,
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
    # TODO: we DO need the contig_seq which can only be read from fasta file
    contigs = {}
    with InputStream(contig_file) as file:
        for rec in Bio.SeqIO.parse(file, 'fasta'):
            contigs[rec.id] = {
                "species_id": species_id,
                "contig_id": rec.id,
                "contig_len": len(rec.seq),
                "contig_seq": str(rec.seq),
            }
    return contigs


def design_chunks(species_ids_of_interest, contigs_files, chunk_size):
    """ Each chunk is indexed by species_id, chunk_id """

    global sample
    global semaphore_for_species
    global species_sliced_snps_path

    semaphore_for_species = dict()
    # design per species list of headerless_sliced_output_path indexed by chunk_id.
    species_sliced_snps_path = defaultdict(list)
    species_sliced_snps_path["input_bamfile"] = sample.get_target_layout("snps_repgenomes_bam")

    # TODO: this part can be done in parallel
    arguments_list = []
    for species_index, species_id in enumerate(species_ids_of_interest):
        # Read in contigs information for one species.
        # TODO: download contigs here. Nope. we need to build one cat-ed bowtie2 database
        contigs = scan_contigs(contigs_files[species_index], species_id)

        chunk_id = 0
        for contig_id in sorted(list(contigs.keys())): # why need to sort?
            contig = contigs[contig_id]
            contig_length = contig["contig_len"]

            if contig_length <= chunk_size:
                headerless_sliced_path = sample.get_target_layout("chunk_pileup", species_id, chunk_id)
                species_sliced_snps_path[species_id].append(headerless_sliced_path)

                # TODO: instead contig as the last argument, just pass the contig_seq.
                slice_args = (species_id, chunk_id, contig_id, 0, contig_length, contig)
                arguments_list.append(slice_args)
                chunk_id += 1
            else:
                number_of_chunks = ceil(contig_length/chunk_size) - 1
                for ni, ci in enumerate(range(0, contig_length, chunk_size)):
                    headerless_sliced_path = sample.get_target_layout("chunk_pileup", species_id, chunk_id)
                    species_sliced_snps_path[species_id].append(headerless_sliced_path)

                    if ni == number_of_chunks:
                        slice_args = (species_id, chunk_id, contig_id, ci, contig_length, contig)
                    else:
                        slice_args = (species_id, chunk_id, contig_id, ci, ci+chunk_size, contig)
                    arguments_list.append(slice_args)
                    chunk_id += 1

        # Submit the merge jobs
        arguments_list.append((species_id, -1))
        species_sliced_snps_path[species_id].append(sample.get_target_layout("snps_pileup", species_id))

        # Create a semaphore with number_of_chunks of elements
        semaphore_for_species[species_id] = multiprocessing.Semaphore(chunk_id)
        for _ in range(chunk_id):
            semaphore_for_species[species_id].acquire()

    return arguments_list


def process_chunk_of_sites(packed_args):

    global semaphore_for_species
    global species_sliced_snps_path

    if packed_args[1] == -1:
        species_id = packed_args[0]
        number_of_chunks = len(species_sliced_snps_path[species_id]) - 1
        for _ in range(number_of_chunks):
            semaphore_for_species[species_id].acquire()
        return merge_chunks_per_species(species_id)

    return compute_pileup_per_chunk(packed_args)


def compute_pileup_per_chunk(packed_args):
    """ actual pileup compute for one chunk """
    global semaphore_for_species
    global species_sliced_snps_path
    global global_args
    args = global_args

    try:
        # [contig_start, contig_end)
        species_id, chunk_id, contig_id, contig_start, contig_end, contig = packed_args
        repgenome_bamfile = species_sliced_snps_path["input_bamfile"]

        headerless_sliced_path = species_sliced_snps_path[species_id][chunk_id]

        zero_rows_allowed = not args.sparse
        current_chunk_size = contig_end - contig_start

        print(packed_args[:4])
        with AlignmentFile(repgenome_bamfile) as bamfile:
            counts = bamfile.count_coverage(contig_id, contig_start, contig_end,
                                            quality_threshold=args.aln_baseq, # min_quality_threshold a base has to reach to be counted.
                                            read_callback=keep_read) # select a call-back to ignore reads when counting
            aligned_reads = bamfile.count(contig_id, contig_start, contig_end)
            mapped_reads = bamfile.count(contig_id, contig_start, contig_end, read_callback=keep_read)

        # aln_stats need to be passed from child process back to parents
        aln_stats = {
            "species_id": species_id,
            "contig_id": contig_id,
            "chunk_length": current_chunk_size,
            "aligned_reads": aligned_reads,
            "mapped_reads": mapped_reads,
            "contig_total_depth": 0,
            "contig_covered_bases": 0
        }

        with OutputStream(headerless_sliced_path) as stream:
            for within_chunk_index in range(0, current_chunk_size):
                depth = sum([counts[nt][within_chunk_index] for nt in range(4)])
                count_a = counts[0][within_chunk_index]
                count_c = counts[1][within_chunk_index]
                count_g = counts[2][within_chunk_index]
                count_t = counts[3][within_chunk_index]

                ref_pos = within_chunk_index + contig_start
                ref_allele = contig["contig_seq"][ref_pos]
                row = (contig_id, ref_pos + 1, ref_allele, depth, count_a, count_c, count_g, count_t)

                aln_stats["contig_total_depth"] += depth
                if depth > 0:
                    aln_stats["contig_covered_bases"] += 1
                if depth > 0 or zero_rows_allowed:
                    stream.write("\t".join(map(format_data, row)) + "\n")
            assert within_chunk_index+contig_start == contig_end-1, f"compute_pileup_per_chunk::index mismatch error for {contig_id}."
        return aln_stats
    finally:
        semaphore_for_species[species_id].release() # no deadlock


def merge_chunks_per_species(species_id):
    """ merge the pileup results from chunks into one file per species """
    global species_sliced_snps_path
    global semaphore_for_species
    global global_args

    tsprint(f"merge_chunks_per_species::{species_id}")
    files_of_chunks = species_sliced_snps_path[species_id][:-1]
    species_snps_pileup_file = species_sliced_snps_path[species_id][-1]

    with OutputStream(species_snps_pileup_file) as stream:
        stream.write('\t'.join(snps_pileup_schema.keys()) + '\n')
    cat_files(files_of_chunks, species_snps_pileup_file, 20)

    if not global_args.debug:
        tsprint(f"Deleting temporary sliced pileup files for {species_id}.")
        for s_file in files_of_chunks:
            command(f"rm -rf {s_file}", quiet=True)
    # return a status flag
    # the path should be computable somewhere else
    return True


def write_species_pileup_summary(chunks_pileup_summary, outfile):
    """ Collect species pileup aln stats from all chunks and write to file """

    species_pileup_summary = defaultdict(dict)
    prev_species_id = None

    for record in chunks_pileup_summary:
        if record is True:
            continue

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
        assert curr_species_pileup["fraction_covered"] == 0.0, f"compute_species_pileup_summary error for {species_id}"

        if prev_species_id != species_id:
            if prev_species_id is not None:
                previous_species_pileup = species_pileup_summary.get(prev_species_id)
                if previous_species_pileup["genome_length"] > 0:
                    previous_species_pileup["fraction_covered"] = previous_species_pileup["covered_bases"] / previous_species_pileup["genome_length"]
                if previous_species_pileup["covered_bases"] > 0:
                    previous_species_pileup["mean_coverage"] = previous_species_pileup["total_depth"] / previous_species_pileup["covered_bases"]
            prev_species_id = species_id

    # Secondary level compute: need to loop over all chunks of sites to compute fraction_covered
    if curr_species_pileup["genome_length"] > 0:
        curr_species_pileup["fraction_covered"] = curr_species_pileup["covered_bases"] / curr_species_pileup["genome_length"]
    if curr_species_pileup["covered_bases"] > 0:
        curr_species_pileup["mean_coverage"] = curr_species_pileup["total_depth"] / curr_species_pileup["covered_bases"]

    with OutputStream(outfile) as stream:
        stream.write("\t".join(snps_profile_schema.keys()) + "\n")
        for record in species_pileup_summary.values():
            stream.write("\t".join(map(format_data, record.values())) + "\n")


def midas_run_snps(args):

    try:
        global sample
        sample = Sample(args.sample_name, args.midas_outdir, "snps")
        sample.create_dirs(["outdir", "tempdir", "dbsdir"], args.debug)

        global global_args
        global_args = args

        species_list = args.species_list.split(",") if args.species_list else []
        if args.local_bowtie2_indexes:
            bt2_db_dir = os.path.dirname(args.local_bowtie2_indexes)
            bt2_db_name = os.path.basename(args.local_bowtie2_indexes)
            assert bowtie2_index_exists(bt2_db_dir, bt2_db_name), f"Provided {bt2_db_dir}/{bt2_db_name} don't exist."
            assert (args.species_profile_path and os.path.exists(args.species_profile_path)), f"Need to provide valid species_profile_path."

            bt2_species = []
            with InputStream(args.species_profile_path) as stream:
                for species_id in select_from_tsv(stream, ["species_id"]):
                    bt2_species.append(species_id[0])
            species_list = list(set(species_list) & set(bt2_species))
            # The index was only for the purpose of same bowtie2 index. but the species_ids_of_interest per sample
            # can still be based on sample itself.
            # We also don't want too many empty species in our parsing stageself. need to fix the species_prevalence.tsv with SampleID??
            # TODO: should we also provide symlink?
        else:
            bt2_db_dir = sample.get_target_layout("dbsdir")
            bt2_db_name = "repgenomes"

        species_ids_of_interest = sample.select_species(args.genome_coverage, species_list)

        # Download per-species UHGG file into temporary dbs directory
        local_toc = download_reference(outputs.genomes, sample.get_target_layout("dbsdir"))
        sample.create_species_subdirs(species_ids_of_interest, "dbstemp", args.debug)
        contigs_files = UHGG(local_toc).fetch_files(species_ids_of_interest, sample.get_target_layout("dbs_tempdir"), filetype="contigs")

        # Build one bowtie database for species in the restricted species profile
        if bowtie2_index_exists(bt2_db_dir, bt2_db_name):
            build_bowtie2_db(bt2_db_dir, bt2_db_name, contigs_files)
        # Perhaps avoid this giant conglomerated file, fetching instead submaps for each species.
        # TODO: Also colocate/cache/download in master for multiple slave subcommand invocations

        # Map reads to the existing bowtie2 indexes
        sample.create_species_subdirs(species_ids_of_interest, "temp", args.debug)
        repgenome_bamfile = sample.get_target_layout("snps_repgenomes_bam")
        bowtie2_align(bt2_db_dir, bt2_db_name, repgenome_bamfile, args)
        samtools_index(repgenome_bamfile, args.debug)

        # Use mpileup to call SNPs
        arguments_list = design_chunks(species_ids_of_interest, contigs_files, args.chunk_size)
        chunks_pileup_summary = multiprocessing_map(process_chunk_of_sites, arguments_list, num_physical_cores)

        write_species_pileup_summary(chunks_pileup_summary, sample.get_target_layout("snps_summary"))
    except:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error. Specify --debug flag to keep.")
            sample.remove_dirs(["outdir", "tempdir", "dbsdir", "dbs_tempdir"])
        raise


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_run_snps(args)
