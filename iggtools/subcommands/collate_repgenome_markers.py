import os
import sys
from collections import defaultdict
from multiprocessing import Semaphore
from iggtools.common.argparser import add_subcommand, SUPPRESS
from iggtools.common.utils import tsprint, InputStream, OutputStream, download_reference, select_from_tsv, retry, command, multithreading_map, multiprocessing_map, cat_files, find_files, upload, num_physical_cores, split, upload_star, pythonpath
from iggtools.models.uhgg import UHGG, MIDAS_IGGDB, MARKER_FILE_EXTS, get_uhgg_layout
from iggtools.params.schemas import MARKER_INFO_SCHEMA, PAN_GENE_INFO_SCHEMA
from iggtools.subcommands.build_pangenome import decode_species_arg
from iggtools.params import outputs


CONCURRENT_SPECIES_BUILDS = Semaphore(32)


def pangenome_file(species_id, component):
    # s3://microbiome-igg/2.0/pangenomes/GUT_GENOMEDDDDDD/{genes.ffn, centroids.ffn, gene_info.txt}
    return f"{outputs.pangenomes}/{species_id}/{component}"


def species_mc_file(species_id, component):
    #.txt.lz4
    return f"{outputs.marker_centroids}/{species_id}.{component}"


def genome_marker_mapfile(species_id, genome_id, component):
    # markers.map.lz4
    return f"{outputs.marker_genes}/temp/{species_id}/{genome_id}/{genome_id}.{component}"


@retry
def find_files_with_retry(f):
    return find_files(f)


def collate_repgenome_markers(args):
    """ Collate marker genes of repgenomes into phyeco.fa and phyeco.map """
    midas_iggdb = MIDAS_IGGDB(args.midas_iggdb)
    species = midas_iggdb.uhgg.species

    collate_log_remote = midas_iggdb.get_target_layout("marker_collate_log", remote=True)
    msg = f"Collating marker genes sequences."
    if find_files_with_retry(collate_log_remote):
        if not args.force:
            tsprint(f"Destination {collate_log_remote} already exists.  Specify --force to overwrite.")
            return
        msg = msg.replace(msg.split(" ")[0], "Re-" + msg.split(" ")[0])
    tsprint(msg)

    collate_log = midas_iggdb.get_target_layout("marker_collate_log", remote=False)
    collate_subdir = os.path.dirname(collate_log)
    if not os.path.isdir(collate_subdir):
        command(f"mkdir -p {collate_subdir}")
    with open(collate_log, "w") as slog:
        slog.write(msg + "\n")

    # Fetch marker genes fasta-file and map-file from s3
    marker_genes_fasta = midas_iggdb.fetch_files("marker_genes_fa", species.keys())
    marker_genes_maps = midas_iggdb.fetch_files("marker_genes_map", species.keys())

    # Collate to phyeco.fa and phyeco.map
    collated_genes_fa = midas_iggdb.get_target_layout("marker_db", remote=False, component="fa")
    for marker_fa_files in split(marker_genes_fasta.values(), 20):
        command("cat " + " ".join(marker_fa_files) + f" >> {collated_genes_fa}")

    collaged_genes_map = midas_iggdb.get_target_layout("marker_db", remote=False, component="map")
    for marker_map_files in split(marker_genes_maps.values(), 20):
        command("cat " + " ".join(marker_map_files) + f" >> {collaged_genes_map}")

    # Build hs-blastn index for the collated phyeco sequences
    cmd_index = f"hs-blastn index {collated_genes_fa} &>> {collate_log}"
    with open(f"{collate_log}", "a") as slog:
        slog.write(cmd_index + "\n")
    command(cmd_index)

    # Upload generated fasta and index files
    upload_tasks = []
    for ext in MARKER_FILE_EXTS:
        local_file = midas_iggdb.get_target_layout("marker_db", remote=False, component=ext)
        s3_file = midas_iggdb.get_target_layout("marker_db", remote=True, component=ext)
        upload_tasks.append((local_file, s3_file))
    multithreading_map(upload_star, upload_tasks)

    # Upload the log file in the last
    upload(collate_log, collate_log_remote, check=False)


def parse_marker_by_genome(packed_ids):
    species_id, genome_id = packed_ids
    marker_by_genome = dict()
    with InputStream(genome_marker_mapfile(species_id, genome_id, "markers.map.lz4")) as stream:
        for gene_id, marker_id in select_from_tsv(stream, ["gene_id", "marker_id"], schema=MARKER_INFO_SCHEMA):
            marker_by_genome[gene_id] = marker_id
    return marker_by_genome


def build_marker_centroids_mapping(args):
    if args.zzz_slave_toc:
        build_mc_slave(args)
    else:
        build_mc_master(args)


def build_mc_master(args):

    # Fetch table of contents from s3.
    # This will be read separately by each species build subcommand, so we make a local copy.
    local_toc = download_reference(outputs.genomes)

    db = UHGG(local_toc)
    species = db.species

    #midas_iggdb = MIDAS_IGGDB(args.midas_iggdb)
    #species = midas_iggdb.uhgg.species

    def species_work(species_id):
        assert species_id in species, f"Species {species_id} is not in the database."
        species_genomes = species[species_id]

        dest_file = species_mc_file(species_id, "txt.lz4")
        msg = f"Building marker genes to centroids for species {species_id} with {len(species_genomes)} total genomes."
        if find_files_with_retry(dest_file):
            if not args.force:
                tsprint(f"Destination {dest_file} for species {species_id} pangenome already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Building", "Rebuilding")

        with CONCURRENT_SPECIES_BUILDS:
            tsprint(msg)
            slave_log = "mc_build.log"
            slave_subdir = f"marker_centroids_99/{species_id}"
            if not args.debug:
                command(f"rm -rf {slave_subdir}")
            if not os.path.isdir(slave_subdir):
                command(f"mkdir -p {slave_subdir}")
            # Recurisve call via subcommand.  Use subdir, redirect logs.
            slave_cmd = f"cd {slave_subdir}; PYTHONPATH={pythonpath()} {sys.executable} -m iggtools collate_repgenome_markers --map_marker_to_centroids -s {species_id} --zzz_slave_mode --zzz_slave_toc {os.path.abspath(local_toc)} {'--debug' if args.debug else ''} &>> {slave_log}"

            with open(f"{slave_subdir}/{slave_log}", "w") as slog:
                slog.write(msg + "\n")
                slog.write(slave_cmd + "\n")
            try:
                command(slave_cmd)
            finally:
                # Cleanup should not raise exceptions of its own, so as not to interfere with any
                # prior exceptions that may be more informative.  Hence check=False.
                upload(f"{slave_subdir}/{slave_log}", species_mc_file(species_id, "log.lz4"), check=False)
                if not args.debug:
                    command(f"rm -rf {slave_subdir}", check=False)

    # Check for destination presence in s3 with up to 10-way concurrency.
    # If destination is absent, commence build with up to 3-way concurrency as constrained by CONCURRENT_SPECIES_BUILDS.
    species_id_list = decode_species_arg(args, species)
    multithreading_map(species_work, species_id_list, 36)


def build_mc_slave(args):
    violation = "Please do not call build_mc_slave directly.  Violation"
    assert args.zzz_slave_mode, f"{violation}:  Missing --zzz_slave_mode arg."
    assert os.path.isfile(args.zzz_slave_toc), f"{violation}: File does not exist: {args.zzz_slave_toc}"

    db = UHGG(args.zzz_slave_toc)
    species = db.species
    species_id = args.species

    assert species_id in species, f"{violation}: Species {species_id} is not in the database."

    species_genomes = species[species_id]
    species_genomes_ids = species_genomes.keys()

    dest_file = species_mc_file(species_id, "txt.lz4")
    local_file = f"{species_id}.txt"

    command(f"aws s3 rm --recursive {dest_file}")

    # Read in all genome's phyeco.mapfile
    list_marker_dicts = multiprocessing_map(parse_marker_by_genome, ((species_id, genome_id) for genome_id in species_genomes_ids))
    # flatten the list of dict
    markers = {k: v for md in list_marker_dicts for k, v, in md.items()}

    # Read in species/genomes/mapfile
    gene_info_file = pangenome_file(species_id, "gene_info.txt.lz4")

    # Filter gene_info.txt by marker_ids
    # Option: can also write to dest_file directly with Outputstream
    with OutputStream(local_file) as ostream:
        ostream.write("\t".join(["marker_id"] + list(PAN_GENE_INFO_SCHEMA.keys())) + "\n")
        with InputStream(gene_info_file) as stream:
            # Loop over the ALl the genes, locate those that are marker genes, and write the
            for row in select_from_tsv(stream, selected_columns=PAN_GENE_INFO_SCHEMA, result_structure=dict):
                if row["gene_id"] in markers.keys():
                    ostream.write("\t".join([markers[row["gene_id"]]] + list(row.values())) + "\n")

    # Update to s3
    upload(local_file, dest_file, check=False)


def register_args(main_func):
    subparser = add_subcommand('collate_repgenome_markers', main_func, help='collate marker genes for repgresentative genomes')
    subparser.add_argument('--midas_iggdb',
                           dest='midas_iggdb',
                           type=str,
                           required=False,
                           metavar="CHAR",
                           help=f"local MIDAS DB which mirrors the s3 IGG db")
    subparser.add_argument('-s',
                           '--species',
                           dest='species',
                           required=False,
                           help="species[,species...] whose pangenome(s) to build;  alternatively, species slice in format idx:modulus, e.g. 1:30, meaning build species whose ids are 1 mod 30; or, the special keyword 'all' meaning all species")
    subparser.add_argument('--collate_repgenome_markers',
                           action='store_true',
                           default=False,
                           help=f"collate marker genes and mapfile for the repgenomes")
    subparser.add_argument('--map_marker_to_centroids',
                           action='store_true',
                           default=False,
                           help=f"Identify the centroids genes for rep marker genes")
    subparser.add_argument('--zzz_slave_toc',
                           dest='zzz_slave_toc',
                           required=False,
                           help=SUPPRESS) # "reserved to pass table of contents from master to slave"
    return main_func


@register_args
def main(args):
    tsprint(f"Executing iggtools subcommand {args.subcommand} with args {vars(args)}.")
    if args.collate_repgenome_markers:
        collate_repgenome_markers(args)
    if args.map_marker_to_centroids:
        build_marker_centroids_mapping(args)
