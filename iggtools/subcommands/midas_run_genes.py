#!/usr/bin/env python3
import json
import os
import multiprocessing

from collections import defaultdict
import numpy as np
import Bio.SeqIO
from pysam import AlignmentFile  # pylint: disable=no-name-in-module

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, InputStream, OutputStream, select_from_tsv, command, multiprocessing_map, multithreading_map, num_physical_cores, cat_files
from iggtools.common.bowtie2 import build_bowtie2_db, bowtie2_align, samtools_index, bowtie2_index_exists, _keep_read
from iggtools.models.uhgg import MIDAS_IGGDB
from iggtools.params.schemas import genes_summary_schema, genes_coverage_schema, format_data
from iggtools.models.sample import Sample


DEFAULT_ALN_COV = 0.75
DEFAULT_MARKER_DEPTH = 3.0
DEFAULT_ALN_MAPID = 94.0
DEFAULT_ALN_READQ = 20
DEFAULT_ALN_MAPQ = 0
DEFAULT_CHUNK_SIZE = 5000


def register_args(main_func):
    subparser = add_subcommand('midas_run_genes', main_func, help='Metagenomic pan-genome profiling')

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

    subparser.add_argument('--marker_depth',
                           type=float,
                           dest='marker_depth',
                           metavar='FLOAT',
                           default=DEFAULT_MARKER_DEPTH,
                           help=f"Include species with >X marker coverage ({DEFAULT_MARKER_DEPTH})")
    subparser.add_argument('--species_list',
                           dest='species_list',
                           type=str,
                           metavar="CHAR",
                           help=f"Comma separated list of species ids")

    #  Alignment flags (bowtie, or postprocessing)
    subparser.add_argument('--aln_speed',
                           type=str,
                           dest='aln_speed',
                           default='very-sensitive',
                           choices=['very-fast', 'fast', 'sensitive', 'very-sensitive'],
                           help='Alignment speed/sensitivity (very-sensitive).  If aln_mode is local (default) this automatically issues the corresponding very-sensitive-local, etc flag to bowtie2.')
    subparser.add_argument('--aln_mode',
                           type=str,
                           dest='aln_mode',
                           default='local',
                           choices=['local', 'global'],
                           help='Global/local read alignment (Default local, corresponds to the bowtie2 --local; global corresponds to the bowtie2 --end-to-end).')
    subparser.add_argument('--aln_interleaved',
                           action='store_true',
                           default=False,
                           help='FASTA/FASTQ file in -1 are paired and contain forward AND reverse reads')

    # Coverage flags
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
                           help=f"Discard reads with DEFAULT_ALN_MAPQ < MAPQ. ({DEFAULT_ALN_MAPQ})")
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
    return main_func


def keep_read(aln):
    global global_args
    args = global_args
    return _keep_read(aln, args.aln_mapid, args.aln_readq, args.aln_mapq, args.aln_cov)


def design_chunks(species_ids_of_interest, centroids_files, chunk_size):
    """ Chunks_of_genes and each chunk is indexed by (species_id, chunk_id) """

    tsprint("CZ::design_chunks::start species_counts: %s" % (len(species_ids_of_interest)))

    global sample
    global semaphore_for_species
    global species_sliced_genes_path
    global species_gene_length

    # Read-only global variables
    semaphore_for_species = dict()
    # For each species, list-of-headerless_sliced_output_path is indexed by chunk_id
    species_sliced_genes_path = defaultdict(list)
    species_sliced_genes_path["input_bamfile"] = sample.get_target_layout("genes_pangenomes_bam")

    species_gene_length = defaultdict(dict)

    arguments_list = []
    for species_id in species_ids_of_interest:
        tsprint(f"  CZ::design_chunks::{species_id}::start for loop ")

        gene_count = 0
        chunk_id = 0
        curr_chunk_genes_dict = defaultdict()
        with InputStream(centroids_files[species_id]) as file:
            # TODO: we should generate the centroids_info.txt
            # while the gene_length should be merged with genes_info for next round of database build
            for centroid in Bio.SeqIO.parse(file, 'fasta'):

                if not chunk_id*chunk_size <= gene_count < (chunk_id+1)*chunk_size:
                    # For each chunk, we need the dict to keep track of the gene_length separately
                    headerless_gene_coverage_path = sample.get_target_layout("chunk_coverage", species_id, chunk_id)
                    arguments_list.append((species_id, chunk_id))
                    species_sliced_genes_path[species_id].append(headerless_gene_coverage_path)
                    species_gene_length[species_id][chunk_id] = curr_chunk_genes_dict

                    chunk_id += 1
                    curr_chunk_genes_dict = defaultdict()

                curr_chunk_genes_dict[centroid.id] = len(centroid.seq)
                gene_count += 1

            headerless_gene_coverage_path = sample.get_target_layout("chunk_coverage", species_id, chunk_id)
            arguments_list.append((species_id, chunk_id))
            species_sliced_genes_path[species_id].append(headerless_gene_coverage_path)
            species_gene_length[species_id][chunk_id] = curr_chunk_genes_dict
            chunk_id += 1
        tsprint(f"  CZ::design_chunks::{species_id}::finish for loop with {chunk_id} chunks")

        # Submit merge tasks for all chunks per species
        arguments_list.append((species_id, -1))
        species_sliced_genes_path[species_id].append(sample.get_target_layout("genes_coverage", species_id))

        semaphore_for_species[species_id] = multiprocessing.Semaphore(chunk_id)
        for _ in range(chunk_id):
            semaphore_for_species[species_id].acquire()

    tsprint("CZ::design_chunks::finish species_counts: %s" % (len(species_ids_of_interest)))
    return arguments_list


def process_chunk_of_genes(packed_args):
    """ Compute coverage of pangenome for given species_id and write results to disk """

    global semaphore_for_species
    global species_sliced_genes_path

    if packed_args[1] == -1:
        species_id = packed_args[0]
        number_of_chunks = len(species_sliced_genes_path[species_id]) - 1
        tsprint(f"  CZ::process_chunk_of_genes::{species_id}::wait merge_chunks_per_species")
        for _ in range(number_of_chunks):
            semaphore_for_species[species_id].acquire()
        tsprint(f"  CZ::process_chunk_of_genes::{species_id}::start merge_chunks_per_species")
        ret = merge_chunks_per_species(species_id)
        tsprint(f"  CZ::process_chunk_of_genes::{species_id}::finish merge_chunks_per_species")
        return ret

    species_id = packed_args[0]
    chunk_id = packed_args[1]
    tsprint(f"  CZ::process_chunk_of_genes::{species_id}-{chunk_id}::start compute_coverage_per_chunk")
    ret = compute_coverage_per_chunk(packed_args)
    tsprint(f"  CZ::process_chunk_of_genes::{species_id}-{chunk_id}::finish compute_coverage_per_chunk")
    return ret


def compute_coverage_per_chunk(packed_args):
    """ Count number of bp mapped to each pan-gene. """

    global semaphore_for_species
    species_id = packed_args[0]
    chunk_id = packed_args[1]
    tsprint(f"    CZ::compute_coverage_per_chunk::{species_id}-{chunk_id}::start")

    try:
        global species_sliced_genes_path
        global species_gene_length

        species_id, chunk_id = packed_args
        pangenome_bamfile = species_sliced_genes_path["input_bamfile"]

        headerless_gene_coverage_path = species_sliced_genes_path[species_id][chunk_id]
        gene_length_dict = species_gene_length[species_id][chunk_id]

        # Statistics needed to be accmulated within each chunk
        chunk_of_gene_ids = sorted(list(gene_length_dict.keys()))
        chunk_genome_size = 0
        chunk_num_covered_genes = 0
        chunk_nz_gene_depth = 0
        chunk_aligned_reads = 0
        chunk_mapped_reads = 0

        with OutputStream(headerless_gene_coverage_path) as stream:
            with AlignmentFile(pangenome_bamfile) as bamfile:
                for gene_id in chunk_of_gene_ids:
                    # Basic compute unit for each gene
                    gene_length = gene_length_dict[gene_id]
                    aligned_reads = bamfile.count(gene_id)
                    mapped_reads = bamfile.count(gene_id, read_callback=keep_read)
                    gene_depth = sum((len(aln.query_alignment_sequence) / gene_length for aln in bamfile.fetch(gene_id)))

                    chunk_genome_size += 1
                    if gene_depth == 0: # Sparse by default.
                        continue

                    chunk_num_covered_genes += 1
                    chunk_nz_gene_depth += gene_depth
                    chunk_aligned_reads += aligned_reads
                    chunk_mapped_reads += mapped_reads

                    vals = [gene_id, gene_length, aligned_reads, mapped_reads, gene_depth, 0.0]
                    stream.write("\t".join(map(format_data, vals)) + "\n")

        current_chunk_size = len(chunk_of_gene_ids)
        tsprint(f"    CZ::compute_coverage_per_chunk::{species_id}-{chunk_id}::finish with nz-genes {chunk_num_covered_genes} out of chunk_size {current_chunk_size}")

        return {
            "species_id": species_id,
            "chunk_id": chunk_id,
            "chunk_genome_size": chunk_genome_size,
            "chunk_num_covered_genes": chunk_num_covered_genes,
            "chunk_nz_gene_depth": chunk_nz_gene_depth,
            "chunk_aligned_reads": chunk_aligned_reads,
            "chunk_mapped_reads": chunk_mapped_reads
        }

    finally:
        semaphore_for_species[species_id].release()


def merge_chunks_per_species(species_id):
    """ Compute coverage of pangenome for given species_id and write results to disk """

    tsprint(f"    CZ::merge_chunks_per_species::{species_id}::start")

    global semaphore_for_species
    global species_sliced_genes_path
    global global_args
    global marker_centroids_files # we can deal with marker_centroids_files global or not later

    all_chunks = species_sliced_genes_path[species_id][:-1]
    species_gene_coverage_path = species_sliced_genes_path[species_id][-1]

    # Extract marker gene reads coverage from all the non-zero covered genes.
    tsprint(f"      CZ::get_marker_coverage_from_chunk::{species_id}::start")
    with InputStream(marker_centroids_files[species_id]) as stream:
        centroids_of_marker = dict(select_from_tsv(stream, selected_columns=["marker_id", "centroid_99"]))
    mc_genes = list(centroids_of_marker.values())
    # TODO note: when centroid_70, multiple marker genes may correspond to one centroid_70 (T OR F)?
    # then that centroin_70 would be single copy anymore.

    # Construct awk command for match marker_genes
    pat_str = " || ".join([f"$1==\"{g}\"" for g in mc_genes])
    awk_command = "awk \'%s\'" % pat_str

    marker_genes_depth = dict(zip(mc_genes, [0.0]*len(mc_genes)))
    args = []
    for chunk_file in all_chunks:
        args.append((chunk_file, awk_command, marker_genes_depth))
    multithreading_map(get_marker_coverage_from_chunk, args, 4)
    tsprint(f"      CZ::get_marker_coverage_from_chunk::{species_id}::finish")

    # Overwrite the chunk_gene_coverage file with updated copy_number
    tsprint(f"      CZ::rewrite_chunk_coverage_file::{species_id}::start")
    median_marker_depth = np.median(list(marker_genes_depth.values()))
    if median_marker_depth > 0:
        args = []
        for chunk_file in all_chunks:
            args.append((chunk_file, median_marker_depth))
        multithreading_map(rewrite_chunk_coverage_file, args, 4)
    tsprint(f"      CZ::rewrite_chunk_coverage_file::{species_id}::finish")

    # Merge chunks' results to files genes_coverage
    with OutputStream(species_gene_coverage_path) as stream:
        stream.write('\t'.join(genes_coverage_schema.keys()) + '\n')
    cat_files(all_chunks, species_gene_coverage_path, 20)
    # TODO: multithreading try to write to same file?

    if not global_args.debug:
        tsprint(f"Deleting temporary sliced coverage files for {species_id}.")
        for s_file in all_chunks:
            command(f"rm -rf {s_file}", quiet=True)

    tsprint(f"    CZ::merge_chunks_per_species::{species_id}::finish")
    return {"species_id": species_id, "chunk_id": -1, "median_marker_depth": median_marker_depth}


def get_marker_coverage_from_chunk(my_args):

    chunk_file, awk_command, marker_genes_depth = my_args
    with InputStream(chunk_file, awk_command) as stream:
        for row in select_from_tsv(stream, schema=genes_coverage_schema, result_structure=dict):
            marker_genes_depth[row["gene_id"]] += row["total_depth"]


def rewrite_chunk_coverage_file(my_args):

    chunk_coverage_path, median_marker_depth = my_args

    c_copies = list(genes_coverage_schema.keys()).index("copy_number")
    c_depth = list(genes_coverage_schema.keys()).index("total_depth")

    add_cn_to_write = []
    with InputStream(chunk_coverage_path) as stream:
        for line in stream:
            vals = line.rstrip("\n").split("\t")
            # Infer gene copy counts
            vals[c_copies] = float(vals[c_depth]) / median_marker_depth
            add_cn_to_write.append(vals)

    with OutputStream(chunk_coverage_path) as stream:
        for vals in add_cn_to_write:
            stream.write("\t".join(map(format_data, vals)) + "\n")


def write_species_coverage_summary(chunks_gene_coverage, species_genes_coverage_path):

    tsprint(f"  CZ::write_species_coverage_summary::start")
    species_coverage_summary = defaultdict(dict)

    for record in chunks_gene_coverage:
        # for the merge task, we return the marker genes coverage

        species_id = record["species_id"]

        if record["chunk_id"] == -1:
            species_coverage_summary[species_id]["median_marker_depth"] = record["median_marker_depth"]
            continue

        if species_id not in species_coverage_summary:
            species_coverage_summary[species_id] = {
                "species_id": species_id,
                "pangenome_size": 0,
                "num_covered_genes": 0,
                "total_nz_gene_depth": 0.0,
                "aligned_reads": 0,
                "mapped_reads": 0,
                "median_marker_depth": 0
            }

        species_coverage_summary[species_id]["pangenome_size"] += record["chunk_genome_size"]
        species_coverage_summary[species_id]["num_covered_genes"] += record["chunk_num_covered_genes"]
        species_coverage_summary[species_id]["total_nz_gene_depth"] += record["chunk_nz_gene_depth"]
        species_coverage_summary[species_id]["aligned_reads"] += record["chunk_aligned_reads"]
        species_coverage_summary[species_id]["mapped_reads"] += record["chunk_mapped_reads"]

    # Need to loop over all the chunks to calculate the average read depths coverage
    with OutputStream(species_genes_coverage_path) as stream:
        stream.write("\t".join(genes_summary_schema.keys()) + "\n")
        for record in species_coverage_summary.values():
            mean_coverage = record["total_nz_gene_depth"] / record["num_covered_genes"]
            fraction_covered = record["num_covered_genes"] / record["pangenome_size"]
            vals = [record["species_id"], record["pangenome_size"], record["num_covered_genes"], \
                    fraction_covered, mean_coverage, record["aligned_reads"], record["mapped_reads"], record["median_marker_depth"]]
            stream.write("\t".join(map(format_data, vals)) + "\n")
    tsprint(f"  CZ::write_species_coverage_summary::finish")


def midas_run_genes(args):

    try:
        global sample
        sample = Sample(args.sample_name, args.midas_outdir, "genes")
        sample.create_dirs(["outdir", "tempdir"], args.debug)

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
        else:
            sample.create_dirs(["bt2_indexes_dir"], args.debug)
            bt2_db_dir = sample.get_target_layout("bt2_indexes_dir")
            bt2_db_name = "pangenomes"

        # Select abundant species present in the sample for SNPs calling
        species_ids_of_interest = species_list if args.marker_depth == -1 else sample.select_species(args.marker_depth, species_list)
        species_counts = len(species_ids_of_interest)
        assert species_counts > 0, f"No (specified) species pass the marker_depth filter, please adjust the marker_depth or species_list"
        tsprint(species_ids_of_interest)

        # Fetch centroids_99 fastas for each species (multiprocessing)
        midas_iggdb = MIDAS_IGGDB(args.midas_iggdb if args.midas_iggdb else sample.get_target_layout("midas_iggdb_dir"), args.num_cores)
        centroids_files = midas_iggdb.fetch_files("centroids", species_ids_of_interest)
        tsprint(centroids_files)

        # Fetch marker's centroids cluster info per species
        global marker_centroids_files
        marker_centroids_files = midas_iggdb.fetch_files("marker_centroids", species_ids_of_interest)
        tsprint(marker_centroids_files)

        # Build Bowtie indexes for species in the restricted species profile
        tsprint(f"CZ::build_bowtie2_indexes::start")
        if not bowtie2_index_exists(bt2_db_dir, bt2_db_name):
            build_bowtie2_db(bt2_db_dir, bt2_db_name, centroids_files, args.num_cores)
        tsprint(f"CZ::build_bowtie2_indexes::finish ({species_counts}) species counts")

        # Align reads to pangenome database
        tsprint(f"CZ::bowtie2_align::start")
        sample.create_species_subdirs(species_ids_of_interest, "temp", args.debug)
        pangenome_bamfile = sample.get_target_layout("genes_pangenomes_bam")
        bowtie2_align(bt2_db_dir, bt2_db_name, pangenome_bamfile, args)
        samtools_index(pangenome_bamfile, args.debug, args.num_cores)
        tsprint(f"CZ::bowtie2_align::finish")

        # Compute coverage of genes in pangenome database
        tsprint(f"CZ::design_chunks::start")
        arguments_list = design_chunks(species_ids_of_interest, centroids_files, args.chunk_size)
        tsprint(f"CZ::design_chunks::finish")

        tsprint(f"CZ::multiprocessing_map::start")
        chunks_gene_coverage = multiprocessing_map(process_chunk_of_genes, arguments_list, args.num_cores)
        tsprint(f"CZ::multiprocessing_map::finish")

        write_species_coverage_summary(chunks_gene_coverage, sample.get_target_layout("genes_summary"))

    except Exception as error:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error.  Specify --debug flag to keep.")
            sample.remove_dirs(["outdir", "tempdir"])
        if not args.prebuilt_bowtie2_indexes:
            sample.remove_dirs(["bt2_indexes_dir"])
        raise error


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_run_genes(args)
