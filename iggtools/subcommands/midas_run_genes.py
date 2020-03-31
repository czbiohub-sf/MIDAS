import json
import os
import multiprocessing

from collections import defaultdict
import numpy as np
import Bio.SeqIO
from pysam import AlignmentFile  # pylint: disable=no-name-in-module

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, InputStream, OutputStream, select_from_tsv, command, multiprocessing_map, download_reference, cat_files, split, num_physical_cores
from iggtools.params import outputs
from iggtools.common.bowtie2 import build_bowtie2_db, bowtie2_align, samtools_index, bowtie2_index_exists
from iggtools.models.uhgg import UHGG, get_uhgg_layout
from iggtools.params.schemas import genes_summary_schema, genes_coverage_schema, MARKER_INFO_SCHEMA, format_data
from iggtools.models.sample import Sample


DEFAULT_ALN_COV = 0.75
DEFAULT_GENOME_COVERAGE = 3.0
DEFAULT_ALN_MAPID = 94.0
DEFAULT_ALN_READQ = 20
DEFAULT_ALN_MAPQ = 0
DEFAULT_CHUNK_SIZE = 5000


def register_args(main_func):
    subparser = add_subcommand('midas_run_genes', main_func, help='metagenomic pan-genome profiling')
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

    subparser.add_argument('--num_cores',
                           dest='num_cores',
                           type=int,
                           metavar="INT",
                           default=num_physical_cores,
                           help=f"Number of physical cores to use ({num_physical_cores})")

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

    subparser.add_argument('--chunk_size',
                           dest='chunk_size',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_CHUNK_SIZE,
                           help=f"Number of genomic sites for the temporary chunk file  ({DEFAULT_CHUNK_SIZE})")

    #  Alignment flags (bowtie, or postprocessing)
    subparser.add_argument('--aln_cov',
                           dest='aln_cov',
                           default=DEFAULT_ALN_COV,
                           type=float,
                           metavar="FLOAT",
                           help=f"Discard reads with alignment coverage < ALN_COV ({DEFAULT_ALN_COV}).  Values between 0-1 accepted.")
    subparser.add_argument('--aln_readq',
                           dest='aln_readq',
                           type=int,
                           metavar="INT",
                           default=DEFAULT_ALN_READQ,
                           help=f"Discard reads with mean quality < READQ ({DEFAULT_ALN_READQ})")
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
                           help='Global/local read alignment (local, corresponds to the bowtie2 --local).  Global corresponds to the bowtie2 default --end-to-end.')
    subparser.add_argument('--aln_interleaved',
                           action='store_true',
                           default=False,
                           help='FASTA/FASTQ file in -1 are paired and contain forward AND reverse reads')
    subparser.add_argument('--aln_sort',
                           action='store_true',
                           default=True,
                           help=f"Sort BAM file.")

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


def marker_to_centroid_mapping(species_id):
    """ Identify which cluster the marker_genes belong to """
    # TODO: pre-build this part for all species in the uhgg db
    # find marker genes from the centroid_99 genes
    # We also have the marker genes for all the genomes phyeco/temp/genome/genome.marker.map/fa

    tsprint(f"  CZ::marker_to_centroid_mapping::{species_id}::start")
    global sample

    # Get the gene_id - marker_id map
    markers = dict()
    awk_command = f"awk \'$1 == \"{species_id}\"\'"
    marker_genes_mapfile = get_uhgg_layout("marker_genes_mapfile")
    with InputStream(marker_genes_mapfile, awk_command) as stream:
        for gene_id, marker_id in select_from_tsv(stream, ["gene_id", "marker_id"], schema=MARKER_INFO_SCHEMA):
            assert marker_id not in markers, f"marker {marker_id} for species {species_id} corresponds to multiple gene_ids."
            markers[gene_id] = marker_id

    # Get the gene_id to centroid_gene_id map
    # TODO This part can also be done during the database build
    gene_info = dict()
    pangenome_file = get_uhgg_layout("pangenome_file", species_id, "gene_info.txt.lz4")
    with InputStream(pangenome_file) as stream:
        for gene_id, centroids_gene_id in select_from_tsv(stream, ["gene_id", "centroid_99"]):
            if gene_id in markers.keys():
                gene_info[gene_id] = centroids_gene_id

    marker_to_centroid_dict = dict()
    for gene_id, marker_id in markers.items():
        centroid_gene_id = gene_info[gene_id]
        marker_to_centroid_dict[marker_id] = centroid_gene_id
        # I don't think one marker_gene should map to multiple centroids genes
        #marker_to_centroid_dict[marker_id].append(centroid_gene_id)

    # Write marker_to_centroid_dict to file
    with open(sample.get_target_layout("marker_genes_mapping", species_id), "a") as stream:
        for k, v in marker_to_centroid_dict.items():
            stream.write("\t".join([k, v]) + "\n")
    tsprint(f"  CZ::marker_to_centroid_mapping::{species_id}::end")


def design_chunks(species_ids_of_interest, centroids_files, chunk_size):
    tsprint("CZ::design_chunks::start num_of_species: %s" % (len(species_ids_of_interest)))
    global sample
    global semaphore_for_species
    global species_sliced_genes_path
    global species_gene_length
    global species_marker_genes

    semaphore_for_species = dict()
    species_sliced_genes_path = defaultdict(list)
    species_gene_length = defaultdict(dict)
    species_marker_genes = defaultdict(dict)

    species_sliced_genes_path["input_bamfile"] = sample.get_target_layout("genes_pangenomes_bam")

    arguments_list = []
    for species_index, species_id in enumerate(species_ids_of_interest):

        tsprint(f"  CZ::design_chunks::{species_id}::start")
        # Get the list of centroids99 genes that contains marker genes in the cluster
        with InputStream(sample.get_target_layout("marker_genes_mapping", species_id)) as stream:
            marker_to_centroid = dict(select_from_tsv(stream, selected_columns=["marker", "centroid"], schema={"marker":str, "centroid":str}))
        c_markers = list(marker_to_centroid.values()) #<= all we care is the list of centroid genes corresponds to marker genes
        species_marker_genes[species_id] = dict(zip(c_markers, [0.0]*len(c_markers)))

        centroid_file = centroids_files[species_index]

        gene_count = 0
        chunk_id = 0
        curr_chunk_genes_dict = defaultdict()

        tsprint(f"  CZ::design_chunks::{species_id}::2. centroid loop start")
        with InputStream(centroid_file) as file:
            # TODO: we should (if not already) have the centroid_99_gene_info.TSV
            # while the gene_length should be merged with genes_info for next round of database build
            # what we need here are: centroid - length
            for centroid in Bio.SeqIO.parse(file, 'fasta'):

                if not chunk_id*chunk_size <= gene_count < (chunk_id+1)*chunk_size:
                    # For each chunk, we need the dict to keep track of the gene_length separately
                    arguments_list.append((species_id, chunk_id))
                    species_sliced_genes_path[species_id].append(sample.get_target_layout("chunk_coverage", species_id, chunk_id))
                    species_gene_length[species_id][chunk_id] = curr_chunk_genes_dict

                    chunk_id += 1
                    curr_chunk_genes_dict = defaultdict()

                curr_chunk_genes_dict[centroid.id] = len(centroid.seq)
                gene_count += 1

            arguments_list.append((species_id, chunk_id))
            species_sliced_genes_path[species_id].append(sample.get_target_layout("chunk_coverage", species_id, chunk_id))
            species_gene_length[species_id][chunk_id] = curr_chunk_genes_dict
            chunk_id += 1

        tsprint(f"  CZ::design_chunks::{species_id}::finish with {chunk_id} chunks")
        # Submit merge tasks for all chunks per species
        arguments_list.append((species_id, -1))
        species_sliced_genes_path[species_id].append(sample.get_target_layout("genes_coverage", species_id))

        semaphore_for_species[species_id] = multiprocessing.Semaphore(chunk_id)
        for _ in range(chunk_id):
            semaphore_for_species[species_id].acquire()

    tsprint("CZ::design_chunks::finish num_of_species: %s" % (len(species_ids_of_interest)))
    return arguments_list


def process_chunk_of_genes(packed_args):
    """ Compute coverage of pangenome for given species_id and write results to disk """
    species_id = packed_args[0]
    tsprint(f"CZ::process_chunk_of_genes::{species_id}::start")

    global semaphore_for_species
    global species_sliced_genes_path

    if packed_args[1] == -1:
        species_id = packed_args[0]
        number_of_chunks = len(species_sliced_genes_path[species_id]) - 1
        tsprint(f"  CZ::process_chunk_of_genes::{species_id}::wait for all chunks to be processed")
        for _ in range(number_of_chunks):
            semaphore_for_species[species_id].acquire()
        tsprint(f"  CZ::process_chunk_of_genes::{species_id}::call merge_chunks_per_species")
        ret = merge_chunks_per_species(species_id)
        tsprint(f"CZ::process_chunk_of_genes::{species_id}::finish")
        return ret

    chunk_id = packed_args[1]
    tsprint(f"  CZ::process_chunk_of_genes::{species_id}-{chunk_id}::start compute_pileup_per_chunk")
    ret = compute_coverage_per_chunk(packed_args)
    tsprint(f"  CZ::process_chunk_of_genes::{species_id}-{chunk_id}::finish compute_pileup_per_chunk")
    return ret


def compute_coverage_per_chunk(packed_args):
    """ Count number of bp mapped to each pan-gene. """

    species_id = packed_args[0]
    chunk_id = packed_args[1]
    tsprint(f"  CZ::compute_coverage_per_chunk::{species_id}-{chunk_id}::start")


    global semaphore_for_species
    global species_sliced_genes_path
    global species_gene_length
    global species_marker_genes

    try:
        species_id, chunk_id = packed_args
        pangenome_bamfile = species_sliced_genes_path["input_bamfile"]

        headerless_gene_coverage_path = species_sliced_genes_path[species_id][chunk_id]
        marker_genes = species_marker_genes[species_id]
        gene_length_dict = species_gene_length[species_id][chunk_id]

        # Statistics needed to be accmulated within each chunk
        chunk_of_gene_ids = sorted(list(gene_length_dict.keys()))
        chunk_genome_size = 0
        chunk_num_covered_genes = 0
        chunk_nz_gene_depth = 0
        chunk_aligned_reads = 0
        chunk_mapped_reads = 0


        tsprint(f"    CZ::compute_coverage_per_chunk::{species_id}-{chunk_id}::1-compute / accumulate stats over genes")
        with OutputStream(headerless_gene_coverage_path) as stream:
            with AlignmentFile(pangenome_bamfile) as bamfile:
                for gene_id in chunk_of_gene_ids:
                    # core, basic compute unit for each gene
                    gene_length = gene_length_dict[gene_id]
                    aligned_reads = bamfile.count(gene_id)
                    mapped_reads = bamfile.count(gene_id, read_callback=keep_read)
                    gene_depth = sum((len(aln.query_alignment_sequence) / gene_length for aln in bamfile.fetch(gene_id)))

                    if gene_id in marker_genes.keys():
                        marker_genes[gene_id] += gene_depth

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
        tsprint(f"    CZ::compute_coverage_per_chunk::{species_id}-{chunk_id}::2-finish with {chunk_num_covered_genes} covered genes out of total {current_chunk_size}")

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
        tsprint(f"  CZ::compute_coverage_per_chunk::{species_id}-{chunk_id}::finish")


def cat_awk_files(files_of_chunks, species_file, median_marker_depth, number_of_chunks=20):
    for temp_files in split(files_of_chunks, number_of_chunks):
        cat_cmd = "cat " + " ".join(temp_files)
        awk_cmd = f"awk -F \"\t\" 'BEGIN {{OFS=FS; ORS=\"\"}}; {{print $1,$2,$3,$4,$5; printf(\"\t%.3f\\n\", $5/{median_marker_depth})}}'"
        command(f"{cat_cmd} | {awk_cmd} >> {species_file}")


def merge_chunks_per_species(species_id):
    """ Compute coverage of pangenome for given species_id and write results to disk """

    global semaphore_for_species
    global species_sliced_genes_path
    global global_args
    global species_marker_genes

    tsprint(f"  CZ::merge_chunks_per_species::{species_id}::start")
    all_chunks = species_sliced_genes_path[species_id][:-1]
    species_gene_coverage_path = species_sliced_genes_path[species_id][-1]

    marker_genes_depth = species_marker_genes[species_id]
    median_marker_depth = np.median(list(marker_genes_depth.values()))
    print(f"median_marker_depth => {median_marker_depth}")


    # Write current species's gene coverage to file
    with OutputStream(species_gene_coverage_path) as stream:
        stream.write('\t'.join(genes_coverage_schema.keys()) + '\n')
    if median_marker_depth > 0:
        cat_awk_files(all_chunks, species_gene_coverage_path, median_marker_depth, 10)
    else:
        cat_files(all_chunks, species_gene_coverage_path, 20)


    if not global_args.debug:
        tsprint(f"Deleting temporary sliced coverage files for {species_id}.")
        for s_file in all_chunks:
            command(f"rm -rf {s_file}", quiet=True)

    tsprint(f"  CZ::merge_chunks_per_species::{species_id}::finish")
    return True


def write_species_coverage_summary(chunks_gene_coverage, species_genes_coverage_path):

    tsprint(f"  CZ::write_species_coverage_summary::start")
    global species_marker_genes

    species_coverage_summary = defaultdict(dict)
    for record in chunks_gene_coverage:
        if record is True:
            continue

        species_id = record["species_id"]
        if species_id not in species_coverage_summary:
            species_coverage_summary[species_id] = {
                "species_id": species_id,
                "pangenome_size": 0,
                "num_covered_genes": 0,
                "total_nz_gene_depth": 0.0,
                "aligned_reads": 0,
                "mapped_reads": 0,
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
            median_marker_depth = np.median(list(species_marker_genes[species_id].values()))

            vals = [record["species_id"], record["pangenome_size"], record["num_covered_genes"], \
                    fraction_covered, mean_coverage, record["aligned_reads"], record["mapped_reads"], median_marker_depth]
            stream.write("\t".join(map(format_data, vals)) + "\n")
    tsprint(f"  CZ::write_species_coverage_summary::finish")


def midas_run_genes(args):
    try:
        global sample
        tsprint(f"CZ::midas_run_genes::Sample::start")
        sample = Sample(args.sample_name, args.midas_outdir, "genes")
        sample.create_dirs(["outdir", "tempdir", "dbsdir"], args.debug)
        tsprint(f"CZ::midas_run_genes::Sample::start")

        global global_args
        global_args = args

        species_list = args.species_list.split(",") if args.species_list else []
        if args.local_bowtie2_indexes:
            # Already-built bowtie2 indexes
            bt2_db_dir = os.path.dirname(args.local_bowtie2_indexes)
            bt2_db_name = os.path.basename(args.local_bowtie2_indexes)
            assert bowtie2_index_exists(bt2_db_dir, bt2_db_name), f"Provided {bt2_db_dir}/{bt2_db_name} don't exist."
            assert (args.species_profile_path and os.path.exists(args.species_profile_path)), f"Need to provide valid species_profile_path."

            # Update species_list: either particular species of interest or species in the bowtie2 indexes
            bt2_species = []
            with InputStream(args.species_profile_path) as stream:
                for species_id in select_from_tsv(stream, ["species_id"]):
                    bt2_species.append(species_id[0])
            species_list = list(set(species_list) & set(bt2_species)) if species_list else bt2_species

        else:
            sample.create_dirs(["bt2_indexes_dir"], args.debug)
            bt2_db_dir = sample.get_target_layout("bt2_indexes_dir")
            bt2_db_name = "pangenomes"

        tsprint(f"CZ::select_species::start")
        species_ids_of_interest = sample.select_species(args.genome_coverage, species_list)
        species_counts = len(species_ids_of_interest)
        tsprint(f"  CZ::number of species to analyze: {species_counts}")
        tsprint(f"CZ::select_species::finish")

        tsprint(f"CZ::download_reference::start")
        # Download per-species UHGG file into temporary dbs directory
        local_toc = download_reference(outputs.genomes, sample.get_target_layout("dbsdir"))
        sample.create_species_subdirs(species_ids_of_interest, "dbs", args.debug)
        centroids_files = UHGG(local_toc).fetch_files(species_ids_of_interest, sample.get_target_layout("dbsdir"), filetype="centroids")
        tsprint(f"CZ::download_reference::finish")

        tsprint(f"CZ::bowtie2_align::start")
        # Build one bowtie database for species in the restricted species profile
        if not bowtie2_index_exists(bt2_db_dir, bt2_db_name):
            build_bowtie2_db(bt2_db_dir, bt2_db_name, centroids_files)
        tsprint(f"CZ::build_bowtie2_indexes::finish ({species_counts}) species counts")


        # Map reads to pan-genes bowtie2 database
        sample.create_species_subdirs(species_ids_of_interest, "temp", args.debug)
        pangenome_bamfile = sample.get_target_layout("genes_pangenomes_bam")
        bowtie2_align(bt2_db_dir, bt2_db_name, pangenome_bamfile, args)
        tsprint(f"CZ::bowtie2_align::finish")

        tsprint(f"CZ::samtools_index::start")
        samtools_index(pangenome_bamfile, args.debug)
        tsprint(f"CZ::samtools_index::finish")

        # Convert marker_genes to centroid_genes for each species
        # TODO: move this part go UHGG
        tsprint("CZ::start marker_to_centroid_mapping start")
        multiprocessing_map(marker_to_centroid_mapping, species_ids_of_interest, args.num_cores)
        tsprint("CZ::start marker_to_centroid_mapping finish")

        tsprint(f"CZ::design_chunks::start")
        arguments_list = design_chunks(species_ids_of_interest, centroids_files, args.chunk_size)
        tsprint(f"CZ::design_chunks::finish")

        tsprint(f"CZ::multiprocessing map::start")
        chunks_gene_coverage = multiprocessing_map(process_chunk_of_genes, arguments_list, args.num_cores)
        tsprint(f"CZ::multiprocessing map::finish")
        write_species_coverage_summary(chunks_gene_coverage, sample.get_target_layout("genes_summary"))
        tsprint(f"CZ::midas_run_genes::finish")
    except:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error.  Specify --debug flag to keep.")
            sample.remove_dirs(["outdir", "tempdir", "dbsdir", "bt2_indexes_dir"])
        raise


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    tsprint(f"CZ::midas_run_genes::start")
    midas_run_genes(args)
