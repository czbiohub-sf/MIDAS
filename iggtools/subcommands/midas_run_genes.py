import json
import os
import multiprocessing

from collections import defaultdict
import numpy as np
import Bio.SeqIO
from pysam import AlignmentFile  # pylint: disable=no-name-in-module

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, command, InputStream, OutputStream, select_from_tsv, multiprocessing_map, download_reference, num_physical_cores
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
DEFAULT_CHUNK_SIZE = 2000

DECIMALS = ".6f"


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


def scan_centroids(centroid_file, species_id):
    """ parse centroid fasta for given species_id"""
    # TODO：part of the information should already be processed during db build
    centroids = {}
    with InputStream(centroid_file) as file:
        for centroid in Bio.SeqIO.parse(file, 'fasta'):
            centroids[centroid.id] = {
                "centroid_gene_id": centroid.id,
                "species_id": species_id,
                "length": len(centroid.seq),
                "depth": 0.0,
                "aligned_reads": 0,
                "mapped_reads": 0,
                "copies": 0.0,
            }
    return centroids


def marker_to_centroid_mapping(species_id):
    """ Build mapping between marker_genes and centroid_genes """
    global sample

    # This need to be cleaned up
    #centroid_file = sample.get_target_layout("centroid_file", species_id)
    #centroids = scan_centroids(centroid_file, species_id)

    # Get the gene_id - marker_id map
    markers = dict()
    awk_command = f"awk \'$1 == \"{species_id}\"\'"
    marker_genes_mapfile = get_uhgg_layout()["marker_genes_mapfile"]
    with InputStream(marker_genes_mapfile, awk_command) as stream:
        for gene_id, marker_id in select_from_tsv(stream, ["gene_id", "marker_id"], schema=MARKER_INFO_SCHEMA):
            assert marker_id not in markers, f"marker {marker_id} for species {species_id} corresponds to multiple gene_ids."
            markers[gene_id] = marker_id

    # Get the gene_id to centroid_gene_id map
    # TODO This part can also be done during the database build
    gene_info = dict()
    pangenome_file = get_uhgg_layout(species_id, "gene_info.txt.lz4")["pangenome_file"]
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
    with OutputStream(sample.get_target_layout("marker_genes_mapping", species_id)) as stream:
        for k, v in marker_to_centroid_dict.items():
            stream.write("\t".join([k, v]) + "\n")


def write_genes_summary(genes_summary_list, genes_summary_path):
    # Collect species-level stats
    # This part is done on a species level
    # slice by species_id
    if False:
        nz_gene_depth = [gd["depth"] for gd in centroids.values() if gd["depth"] > 0]
        num_covered_genes = len(nz_gene_depth) #count # covered genes
        mean_coverage = np.mean(nz_gene_depth) # average gene depth
        pangenome_size = len(centroids)
        fraction_covered = num_covered_genes / pangenome_size
        aligned_reads = sum(g["aligned_reads"] for g in centroids.values())
        mapped_reads = sum(g["mapped_reads"] for g in centroids.values())

        values = [species_id, pangenome_size, num_covered_genes, fraction_covered, \
                  mean_coverage, marker_median_depth, aligned_reads, mapped_reads]

    with OutputStream(genes_summary_path) as stream:
        stream.write("\t".join(genes_summary_schema.keys()) + "\n")
        for record in genes_summary_list:
            stream.write("\t".join(map(format_data, record)) + "\n")


def compute_and_write_chunks_per_species(species_id):
    """ Compute coverage of pangenome for given species_id and write results to disk """

    global semaphore_for_species
    global species_sliced_genes_path

    all_chunks = species_sliced_genes_path[species_id]
    gene_coverage_path = species_sliced_genes_path[species_id][-1][0]
    marker_to_centroid_path = species_sliced_genes_path[species_id][-1][1]

    marker_to_centroid = {}
    with InputStream(marker_to_centroid_path) as stream:
        for line in stream:
            line = line.rstrip("\n").split("\t")
            marker_to_centroid[line[0]] = line[1]
    c_markers = marker_to_centroid.values() #<= all we care is the list of centroid genes corresponds to marker genes

    # first compute the marker genes median depths
    marker_genes_depth = []
    for chunk_id, chunk_of_genes in all_chunks.items():
        for gene_id, gene_val in chunk_of_genes.items():
            if gene_id in c_markers:
                marker_genes_depth.append(gene_val["depth"])
    species_marker_depth = np.median(marker_genes_depth)
    tsprint(f"species_marker_depth => {species_marker_depth}")

    # Write current species's gene coverage to file
    with OutputStream(gene_coverage_path) as stream:
        stream.write('\t'.join(genes_coverage_schema.keys()) + '\n')
        for chunk_id, chunk_of_genes in all_chunks.items():
            for gene_id, gene_dict in chunk_of_genes.items():
                if gene_dict["depth"] == 0: # Sparse by default here.
                    continue
                # second infer gene copy counts
                if species_marker_depth > 0:
                    gene_dict["copies"] = gene_dict["depth"] / species_marker_depth
                print(gene_dict)
                vals = [gene_id, gene_dict["length"], gene_dict["aligned_reads"], gene_dict["mapped_reads"], gene_dict["depth"], gene_dict["copies"]]
                stream.write("\t".join(map(format_data, vals)) + "\n")
    return True


def process_chunk(packed_args):
    """ Compute coverage of pangenome for given species_id and write results to disk """

    global semaphore_for_species
    global species_sliced_genes_path

    if packed_args[1] == -1:
        species_id = packed_args[0]
        number_of_chunks = len(species_sliced_genes_path[species_id])
        for _ in range(number_of_chunks):
            semaphore_for_species[species_id].acquire()

        # All chunk_of_genes coverage information are saved in species_sliced_genes_path[species_id][chunk_id]
        flag = compute_and_write_chunks_per_species(species_id)
        assert flag == True, f"Error for write genes coverage for {species_id}"
        return "worked"

    return compute_chunk_of_genes_coverage(packed_args)


def compute_chunk_of_genes_coverage(packed_args):
    """ Count number of bp mapped to each pan-gene. """

    global semaphore_for_species
    global species_sliced_genes_path

    try:
        species_id, chunk_id, genelist_file = packed_args
        pangenome_bamfile = species_sliced_genes_path["input_bamfile"]

        with InputStream(genelist_file) as stream:
            for line in stream:
                chunk_of_gene_ids = line.rstrip().split("\t")

        genes_acc = species_sliced_genes_path[species_id][chunk_id]
        with AlignmentFile(pangenome_bamfile) as bamfile:
            for gene_id in sorted(chunk_of_gene_ids):
                gene_length = genes_acc[gene_id]["length"]
                genes_acc[gene_id]["aligned_reads"] = bamfile.count(gene_id)
                genes_acc[gene_id]["mapped_reads"] = bamfile.count(gene_id, read_callback=keep_read)
                genes_acc[gene_id]["depth"] = sum((len(aln.query_alignment_sequence) / gene_length for aln in bamfile.fetch(gene_id)))
        return "worked"
    finally:
        semaphore_for_species[species_id].release()


def design_chunks(species_ids_of_interest, chunk_size):
    global sample
    global semaphore_for_species
    global species_sliced_genes_path

    semaphore_for_species = dict()
    species_sliced_genes_path = defaultdict(dict)
    species_sliced_genes_path["input_bamfile"] = sample.get_target_layout("genes_pangenomes_bam")

    arguments_list = []
    for species_id in species_ids_of_interest:
        centroid_file = sample.get_target_layout("centroid_file", species_id)

        gene_count = 0
        chunk_id = 0
        curr_chunk_genes_dict = defaultdict()
        with InputStream(centroid_file) as file:
            # TODO: during database build, this should be solved ...
            for centroid in Bio.SeqIO.parse(file, 'fasta'):
                #if not (gene_count >= chunk_id*chunk_size and gene_count < (chunk_id+1)*chunk_size):
                    #print(f"current chunk {chunk_id}:{gene_count}")
                    #print(f"start new chunk => {chunk_id}")
                if not (gene_count >= chunk_id*chunk_size and gene_count < (chunk_id+1)*chunk_size):
                    chunk_id += 1

                    species_sliced_genes_path[species_id][chunk_id] = curr_chunk_genes_dict
                    gene_list_path = sample.get_target_layout("genes_list", species_id, chunk_id)
                    with OutputStream(gene_list_path) as stream:
                        stream.write("\t".join(curr_chunk_genes_dict.keys()) + "\n")
                    arguments_list.append((species_id, chunk_id, gene_list_path))
                    curr_chunk_genes_dict = defaultdict()

                curr_chunk_genes_dict[centroid.id] = {
                    "centroid_gene_id": centroid.id,
                    "species_id": species_id,
                    "length": len(centroid.seq),
                    "depth": 0.0,
                    "aligned_reads": 0,
                    "mapped_reads": 0,
                    "copies": 0.0,
                }
                gene_count += 1

            species_sliced_genes_path[species_id][chunk_id] = curr_chunk_genes_dict
            gene_list_path = sample.get_target_layout("genes_list", species_id, chunk_id)
            with OutputStream(gene_list_path) as stream:
                stream.write("\t".join(curr_chunk_genes_dict.keys()) + "\n")
            arguments_list.append((species_id, chunk_id, gene_list_path))

        # Submit merge tasks for all chunks per species
        arguments_list.append((species_id, -1))
        species_sliced_genes_path[species_id][-1] = (sample.get_target_layout("genes_coverage", species_id), sample.get_target_layout("marker_genes_mapping", species_id))

        semaphore_for_species[species_id] = multiprocessing.Semaphore(chunk_id)
        for _ in range(chunk_id):
            semaphore_for_species[species_id].acquire()
    return arguments_list


def midas_run_genes(args):

    global sample
    sample = Sample(args.sample_name, args.midas_outdir, "genes")
    sample.create_output_dir(args.debug)

    global global_sample
    global_sample = sample
    global global_args
    global_args = args

    try:
        if args.bt2_db_indexes:
            assert bowtie2_index_exists(bt2_db_dir, bt2_db_name) and os.path.exists(args.species_profile), f"Check the path bt2_db_dir and exists of species_profile_path"
            tsprint("Prebuild bowtie2 index and species_profile exit. Use them")

            bt2_db_dir = os.path.dirname(args.bt2_db_indexes)
            bt2_db_name = os.path.basename(args.bt2_db_indexes)

            species_profile = {}
            with InputStream(args.species_profile_path) as stream:
                for species_id, coverage in select_from_tsv(stream, ["species_id", "coverage"]):
                    species_profile[species_id] = coverage

            species_ids_of_interest = list(species_profile.keys())
            sample.create_species_subdir(species_profile.keys(), args.debug, "dbs")
        else:
            bt2_db_dir = sample.get_target_layout("dbsdir")
            bt2_db_name = "pangenomes"
            bt2_db_temp_dir = sample.get_target_layout("dbs_tempdir")

            if args.species_list:
                species_profile = sample.select_species(args.genome_coverage, args.species_list)
            else:
                species_profile = sample.select_species(args.genome_coverage)

            species_ids_of_interest = list(species_profile.keys())
            sample.create_species_subdir(species_ids_of_interest, args.debug, "dbs")

            # Build one bowtie database for species in the restricted species profile
            local_toc = download_reference(outputs.genomes, bt2_db_dir)
            db = UHGG(local_toc)
            centroids_files = db.fetch_centroids(species_profile.keys(), bt2_db_temp_dir)
            build_bowtie2_db(bt2_db_dir, bt2_db_name, centroids_files)
            # Perhaps avoid this giant conglomerated file, fetching instead submaps for each species.
            # TODO: Also colocate/cache/download in master for multiple slave subcommand invocations.

        species_ids_of_interest = species_ids_of_interest[:2]
        sample.create_species_subdir(species_ids_of_interest, args.debug, "temp")

        # Map reads to pan-genes bowtie2 database
        pangenome_bamfile = sample.get_target_layout("genes_pangenomes_bam")
        bowtie2_align(bt2_db_dir, bt2_db_name, pangenome_bamfile, args)
        samtools_index(pangenome_bamfile, args.debug)

        # Convert marker_genes to centroid_genes for each species
        multiprocessing_map(marker_to_centroid_mapping, species_ids_of_interest, num_physical_cores)

        arguments_list = design_chunks(species_ids_of_interest, args.chunk_size)
        print(arguments_list)

        results = multiprocessing_map(process_chunk, arguments_list, num_physical_cores)
        print(results)
        exit(0)

        # Pangenome coverage compute
        genes_summary_list = []
        for species_index, species_id in enumerate(species_ids_of_interest):
            coverage_path = sample.get_target_layout("genes_coverage", species_id)
            genes_summary_list.append(compute_species_coverage(species_id, centroids_files[species_index], pangenome_bamfile, coverage_path))

        write_genes_summary(genes_summary_list, sample.get_target_layout("genes_summary"))

    except:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error.  Specify --debug flag to keep.")
            sample.remove_output_dir()
        raise


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


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_run_genes(args)
