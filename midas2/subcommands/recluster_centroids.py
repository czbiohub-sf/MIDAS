#!/usr/bin/env python3
import os
import sys
from collections import defaultdict
from multiprocessing import Semaphore

from midas2.common.argparser import add_subcommand
from midas2.common.utils import tsprint, InputStream, command, hashmap, multithreading_map, select_from_tsv, pythonpath, num_physical_cores, copy_star, flatten
from midas2.common.utilities import decode_species_arg
from midas2.models.midasdb import MIDAS_DB
from midas2.params.inputs import MIDASDB_NAMES
from midas2.subcommands.build_pangenome import vsearch, CLUSTERING_PERCENTS, read_uclust_info, write_gene_info, localpath


"""
Input:
    - temp/cdhit/gene_info.tsv
    - temp/cdhit/centroids.99.ffn
    - genes.len
    - genes.ffn
Output:
    - gene_info.txt
    - centroids.xx.ffn
"""


# Up to this many concurrent species builds.
CONCURRENT_SPECIES_BUILDS = Semaphore(1)


def get_gene_info(centroid_info, gene_info_file, percent_id):
    # Parse intermediate gene_info_cdhit.tsv
    with InputStream(gene_info_file) as stream:
        for r in select_from_tsv(stream, selected_columns=['gene_id', 'centroid_99'], schema={'gene_id':str, 'centroid_99':str}):
            centroid_info[r[0]][percent_id] = r[1]


def xref(cluster_files, gene_info_file):
    # Again, let centroid_info[gene][percent_id] be the centroid of the percent_id
    # cluster containing gene; then reclustered to lower percent_id's.
    # Output: gene_info.txt for filtered centroids_99 genereated by cd-hit

    percents = cluster_files.keys()
    max_percent_id = max(percents)
    centroid_info = defaultdict(dict)
    for percent_id, (_, uclust_file) in cluster_files.items():
        if percent_id == max_percent_id:
            get_gene_info(centroid_info, uclust_file, percent_id)
        else:
            read_uclust_info(centroid_info, uclust_file, percent_id)

    # Check for a problem that occurs with improper import of genomes (when contig names clash).
    for g in centroid_info:
        cg = centroid_info[g][max_percent_id]
        ccg = centroid_info[cg][max_percent_id]
        assert cg == ccg, f"The {max_percent_id}-centroid relation should be idempotent, however, {cg} != {ccg}."

    # Infer coarser clusters assignments for all genes by transitivity
    for gc in centroid_info.values():
        gc_recluster = centroid_info[gc[max_percent_id]]
        for percent_id in percents:
            gc[percent_id] = gc_recluster[percent_id]

    # Write to gene_info.txt
    write_gene_info(centroid_info, percents, gene_info_file)


def recluster_centroid(args):
    if args.zzz_worker_mode:
        recluster_centroid_worker(args)
    else:
        recluster_centroid_master(args)


def recluster_centroid_master(args):

    # Fetch table of contents from s3.
    # This will be read separately by each species build subcommand, so we make a local copy.
    midas_db = MIDAS_DB(os.path.abspath(args.midasdb_dir), args.midasdb_name)
    species = midas_db.uhgg.species
    num_threads = args.num_threads

    def species_work(species_id):
        assert species_id in species, f"Species {species_id} is not in the database."
        species_genomes = species[species_id]

        # The species build will upload this file last, after everything else is successfully uploaded.
        # Therefore, if this file exists locally, there is no need to redo the species build.
        dest_file = midas_db.get_target_layout("pangenome_genes_info", False, species_id)
        msg = f"Reclustering centroids for species {species_id} with {len(species_genomes)} total genomes."

        if os.path.exists(dest_file):
            if not args.force:
                tsprint(f"Destination {dest_file} for species {species_id} pangenome already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Building", "Rebuilding")

        with CONCURRENT_SPECIES_BUILDS:
            tsprint(msg)
            last_dest = midas_db.get_target_layout("recluster_log", True, species_id)
            local_dest = midas_db.get_target_layout("recluster_log", False, species_id)
            local_dir = os.path.dirname(local_dest)
            command(f"mkdir -p {local_dir}/temp/cdhit")

            worker_log = os.path.basename(local_dest)
            worker_subdir = os.path.dirname(local_dest) if args.scratch_dir == "." else f"{args.scratch_dir}/buildpan/{species_id}"
            worker_log = f"{worker_subdir}/{worker_log}"
            if not args.debug:
                command(f"rm -rf {worker_subdir}")

            worker_subdir = f"{worker_subdir}/temp"
            if not os.path.isdir(f"{worker_subdir}/cdhit"):
                command(f"mkdir -p {worker_subdir}/cdhit")

            # Recurisve call via subcommand.  Use subdir, redirect logs.
            subcmd_str = f"--zzz_worker_mode -t {num_threads} --midasdb_name {args.midasdb_name} --midasdb_dir {os.path.abspath(args.midasdb_dir)} {'--debug' if args.debug else ''} {'--upload' if args.upload else ''} --scratch_dir {args.scratch_dir}"
            worker_cmd = f"cd {worker_subdir}; PYTHONPATH={pythonpath()} {sys.executable} -m midas2 recluster_centroids -s {species_id} {subcmd_str} &>> {worker_log}"

            with open(f"{worker_log}", "w") as slog:
                slog.write(msg + "\n")
                slog.write(worker_cmd + "\n")
            try:
                command(worker_cmd)
            finally:
                if not args.debug:
                    command(f"rm -rf {worker_subdir}", check=False)
                if args.scratch_dir != ".":
                    command(f"cp -r {worker_log} {local_dest}")

    # Check for destination presence in s3 with up to 8-way concurrency.
    # If destination is absent, commence build with up to 3-way concurrency as constrained by CONCURRENT_SPECIES_BUILDS.
    species_id_list = decode_species_arg(args, species)
    CONCURRENT_RUNS = min(1, int(args.num_threads / 8))
    multithreading_map(species_work, species_id_list, num_threads=CONCURRENT_RUNS)


def recluster_centroid_worker(args):

    violation = "Please do not call generate_gene_info_worker directly.  Violation"
    assert args.zzz_worker_mode, f"{violation}:  Missing --zzz_worker_mode arg."

    species_id = args.species
    midas_db = MIDAS_DB(args.midasdb_dir, args.midasdb_name)

    species = midas_db.uhgg.species
    assert species_id in species, f"{violation}: Species {species_id} is not in the database."

    # start with existing centroids_99 generated by cd-hit
    max_percent, lower_percents = CLUSTERING_PERCENTS[0], CLUSTERING_PERCENTS[1:]
    cluster_files = {max_percent: ("cdhit/centroids.99.ffn", "cdhit/gene_info.txt")}
    assert os.path.isfile("cdhit/centroids.99.ffn"), f"Error: need to run tidy pipeline."
    assert os.path.isfile("cdhit/gene_info.txt"), f"Error: need to run tidy pipeline."

    # reclustering of the max_percent centroids is usually quick, and can proceed in prallel.
    recluster = lambda percent_id: vsearch(percent_id, cluster_files[max_percent][0], num_threads=args.num_threads)
    cluster_files.update(hashmap(recluster, lower_percents))

    # cluster centroids.99 and generate gene_info.txt
    xref(cluster_files, "gene_info.txt")

    if args.scratch_dir != ".":
        copy_tasks = [
            ("cdhit/genes.ffn", localpath(midas_db, species_id, "genes.ffn")),
            ("cdhit/genes.len", localpath(midas_db, species_id, "genes.len")),
            ("gene_info.txt", localpath(midas_db, species_id, "gene_info.txt")),
            (f"cdhit/centroids.{max_percent}.ffn", midas_db.get_target_layout("pangenome_centroids", False, species_id))
        ]
        for src in flatten(cluster_files.values()):
            copy_tasks.append((src, localpath(midas_db, species_id, f"temp/{src}")))
        multithreading_map(copy_star, copy_tasks, args.num_threads)
    # TODO: add upload task on demand


def register_args(main_func):
    subparser = add_subcommand('recluster_centroids', main_func, help='Generate the Generate variety of MIDAS DB related files desired by MIDAS')
    subparser.add_argument('-s',
                           '--species',
                           dest='species',
                           required=False,
                           help="species[,species...] whose pangenome(s) to build;  alternatively, species slice in format idx:modulus, e.g. 1:30, meaning build species whose ids are 1 mod 30; or, the special keyword 'all' meaning all species")
    subparser.add_argument('--midasdb_name',
                           dest='midasdb_name',
                           type=str,
                           default="uhgg",
                           choices=MIDASDB_NAMES,
                           help="MIDAS Database name.")
    subparser.add_argument('--midasdb_dir',
                           dest='midasdb_dir',
                           type=str,
                           default=".",
                           help="Path to local MIDAS Database.")
    subparser.add_argument('-t',
                           '--num_threads',
                           dest='num_threads',
                           type=int,
                           default=num_physical_cores,
                           help="Number of threads")
    subparser.add_argument('--upload',
                           action='store_true',
                           default=False,
                           help="Upload built files to AWS S3.")
    subparser.add_argument('--scratch_dir',
                           dest='scratch_dir',
                           type=str,
                           default=".",
                           help="Absolute path to scratch directory for fast I/O.")
    return main_func


@register_args
def main(args):
    tsprint(f"Executing midas2 subcommand {args.subcommand}.") # with args {vars(args)}.
    recluster_centroid(args)
