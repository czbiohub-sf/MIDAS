#!/usr/bin/env python3
import os
import sys
from collections import defaultdict
import pandas as pd

from midas2.common.argparser import add_subcommand
from midas2.common.utils import tsprint, find_files, OutputStream, command, multithreading_map, pythonpath, num_physical_cores, upload_star, copy_star
from midas2.common.utilities import decode_species_arg, get_contig_length, render_full_stacked_cxx_by_genome, generate_cluster_xx_info, list_cxx_coordinate_to_genome
from midas2.models.midasdb import MIDAS_DB
from midas2.params.inputs import MIDASDB_NAMES
from midas2.subcommands.build_pangenome import localpath


"""
INPUT::FROM build_pangenome + recluster_centroids
    - genes_info.tsv, produced by recluster_centroids.py
    - centroids.ffn, symlink to temp/cdhit/centroids.99.ffn
OUTPUT:
    - augment/{centroid_xx}_info.tsv: centroid_xx, centroid_xx_length, cluster_prevalence, cluster_gene_counts, cluster_marker_id, etc
    - augment/{centroid_xx}_lifted_to_{rep_genome}.tsv
    - augment/{centroid_xx}_stacked_matrix.tsv
"""


def collect_contigs_x_genomes(midas_db, species_id):
    ffns_by_genomes = midas_db.fetch_files("annotation_fna", [species_id], False)
    dict_of_contig_length = defaultdict(dict)
    for genome_id, genomefna in ffns_by_genomes.items():
        dict_of_contig_length[genome_id] = get_contig_length(genomefna)
    return dict_of_contig_length


def write_contig_length(dict_of_contig_length, contig_length_fp):
    with OutputStream(contig_length_fp) as stream:
        stream.write("\t".join(["genome_id", "contig_id", "contig_length"]) + "\n")
        for gid, r in dict_of_contig_length.items():
            for cid, clen in r.items():
                vals = [gid, cid, str(clen)]
                stream.write("\t".join(vals) + "\n")


def write_cxx_info(c99_df, cluster_info_fp):
    c99_df.to_csv(cluster_info_fp, sep='\t', index=False)


def augment_pangenome(args):
    if args.zzz_worker_mode: #"zzz_worker_mode" in args and args.zzz_worker_mode:
        augment_pangenome_worker(args)
    else:
        augment_pangenome_master(args)


def augment_pangenome_master(args):

    midas_db = MIDAS_DB(os.path.abspath(args.midasdb_dir), args.midasdb_name)
    species = midas_db.uhgg.species
    num_threads = args.num_threads

    def species_work(species_id):
        assert species_id in species, f"Species {species_id} is not in the database."
        species_genomes = species[species_id]
        msg = f"Augmenting cluster_xx_info AND contig_length for species {species_id} with {len(species_genomes)} total genomes."

        dest_file = midas_db.get_target_layout("pangenome_contigs_len", True, species_id)
        local_file = midas_db.get_target_layout("pangenome_contigs_len", False, species_id)
        local_dir = os.path.dirname(local_file)
        command(f"mkdir -p {local_dir}/augment")

        if args.upload and find_files(dest_file):
            if not args.force:
                tsprint(f"Destination {dest_file} for species {species_id} cluster_info already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Augmenting", "Reaugmenting")
        if not args.upload and os.path.exists(local_file):
            if not args.force:
                tsprint(f"Destination {local_file} for species {species_id} cluster_info already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Augmenting", "Reaugmenting")

        tsprint(msg)
        last_dest = midas_db.get_target_layout("augment_log", True, species_id)
        local_dest = midas_db.get_target_layout("augment_log", False, species_id)
        local_dir = os.path.dirname(local_dest)
        command(f"mkdir -p {local_dir}/augment")

        worker_log = os.path.basename(local_dest)
        worker_subdir = os.path.dirname(local_dest) if args.scratch_dir == "." else f"{args.scratch_dir}/buildpan/{species_id}"
        worker_log = f"{worker_subdir}/{worker_log}"
        if not args.debug:
            command(f"rm -rf {worker_subdir}")

        worker_subdir = f"{worker_subdir}/augment"
        if not os.path.isdir(worker_subdir):
            command(f"mkdir -p {worker_subdir}")

        # Recurisve call via subcommand.  Use subdir, redirect logs.
        subcmd_str = f"--zzz_worker_mode -t {num_threads} --midasdb_name {args.midasdb_name} --midasdb_dir {os.path.abspath(args.midasdb_dir)} {'--debug' if args.debug else ''} {'--upload' if args.upload else ''} --scratch_dir {args.scratch_dir}"
        worker_cmd = f"cd {worker_subdir}; PYTHONPATH={pythonpath()} {sys.executable} -m midas2 augment_pangenome --species {species_id} {subcmd_str} &>> {worker_log}"

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
            # TODO: handle upload

    species_id_list = decode_species_arg(args, species)
    multithreading_map(species_work, species_id_list, num_threads)


def augment_pangenome_worker(args):

    violation = "Please do not call augment_pangenome_worker directly.  Violation"
    assert args.zzz_worker_mode, f"{violation}:  Missing --zzz_worker_mode arg."

    species_id = args.species
    midas_db = MIDAS_DB(args.midasdb_dir, args.midasdb_name)

    species = midas_db.uhgg.species
    assert species_id in species, f"{violation}: Species {species_id} is not in the database."
    repgenome = midas_db.uhgg.representatives[species_id]

    genes_info_fp = midas_db.fetch_file("pangenome_genes_info", species_id)
    assert os.path.isfile(genes_info_fp), f"Pangenome outputs {genes_info_fp} are missing."
    genes_info = pd.read_csv(genes_info_fp, delimiter="\t")

    # We need this for within-contig-gene-position
    dict_of_contig_length = collect_contigs_x_genomes(midas_db, species_id)
    write_contig_length(dict_of_contig_length, "contigs.len")

    copy_tasks = [
        ("contigs.len", midas_db.get_target_layout("pangenome_contigs_len", False, species_id)),
    ]

    for xx in ["99", "95", "90", "85", "80", "75"]:
        centroid_xx = f"centroid_{xx}"
        cluster_xx_info = generate_cluster_xx_info(genes_info, centroid_xx)
        write_cxx_info(cluster_xx_info, f"clusters_{xx}_info.tsv")
        if args.scratch_dir != ".":
            copy_tasks.append((f"clusters_{xx}_info.tsv", midas_db.get_target_layout("cluster_xx_info", False, species_id, "", xx)))

    centroid_xx_info = render_full_stacked_cxx_by_genome(genes_info, "centroid_75")
    centroid_xx_info.to_csv(f"centroid_75_stacked_matrix.tsv", sep='\t', index=False)

    cxx_lifted_to_qry = list_cxx_coordinate_to_genome(genes_info, "centroid_99", repgenome)
    cxx_lifted_to_qry.to_csv(f"centroid_99_lifted_to_{repgenome}.tsv", sep='\t', index=False)

    if args.scratch_dir != ".":
        copy_tasks.append(("centroid_75_lifted_to_{repgenome}.tsv", localpath(midas_db, species_id, f"augment/centroid_75_lifted_to_{repgenome}.tsv")))
        copy_tasks.append(("centroid_99_stacked_matrix.tsv", localpath(midas_db, species_id, "augment/centroid_99_stacked_matrix.tsv")))

    if args.upload:
        upload_tasks = [
            ("contigs.len", midas_db.get_target_layout("pangenome_contigs_len", True, species_id)),
        ]
        multithreading_map(upload_star, upload_tasks, 2)

    multithreading_map(copy_star, copy_tasks, 2)
    # Remove temp files
    command("rm -f contigs.len", check=False)


def register_args(main_func):
    subparser = add_subcommand('augment_pangenome', main_func, help='Augment pangenome with prevalence and functional annotation')
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
    augment_pangenome(args)
