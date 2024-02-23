#!/usr/bin/env python3
import os
import sys
import json
import pandas as pd
import Bio.SeqIO

from midas2.common.argparser import add_subcommand
from midas2.common.utils import InputStream, tsprint, command, multithreading_map, pythonpath, num_physical_cores, copy_star
from midas2.common.utilities import decode_species_arg
from midas2.models.midasdb import MIDAS_DB
from midas2.params.inputs import MIDASDB_NAMES
from midas2.subcommands.build_pangenome import localpath


"""
We apply a few filters to the centroid_99s for
(1) reducing the bowtie index size without creating reference bias, and
(2) removing potentially across-species contamination, identified as singleton cluster at ANI 75% for genome_prevalence < 0.01 AND genome_counts == 1

Input:
    - clusters_{xx}_info.tsv from enhance_pangenome
Output:
    - pruned/
"""

DEFAULT_PRUNE_CUTOFF = 0.4

def register_args(main_func):
    subparser = add_subcommand('prune_centroids', main_func, help='Prune centroids for specified genomes/species')
    subparser.add_argument('-s',
                           '--species',
                           dest='species',
                           required=False,
                           help="species[,species...] whose centroid(s) to prune;  alternatively, species slice in format idx:modulus, e.g. 1:30, meaning prune species whose ids are 1 mod 30; or, the special keyword 'all' meaning all species")
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
    subparser.add_argument('--scratch_dir',
                           dest='scratch_dir',
                           type=str,
                           default=".",
                           help="Absolute path to scratch directory for fast I/O.")
    subparser.add_argument('--prune_method',
                           dest='prune_method',
                           type=str,
                           default='max',
                           choices=['max', 'median', 'mean'],
                           help="Prune methods: max, median or mean.")
    subparser.add_argument('--prune_cutoff',
                           dest='prune_cutoff',
                           type=float,
                           default=DEFAULT_PRUNE_CUTOFF,
                           help=f"Prune cutoff: for each centroid_95, centroid99 shorter than {DEFAULT_PRUNE_CUTOFF} of the chosen method are pruned for reading mapping.")
    subparser.add_argument('--remove_singleton',
                           action='store_true',
                           default=False,
                           help='Remove c75 clusters with only one gene member in species with more genomes.')
    return main_func


def local_file(midas_db, species_id, xx):
    return midas_db.get_target_layout("pangenome_cluster_xx", False, species_id, "", xx)


def remove_singleton_c75_ids(c75_info):
    c75_info = c75_info[~((c75_info['centroid_75_gene_counts'] == 1) & (c75_info['centroid_75_genome_prevalence'] < 0.1) & (c75_info['centroid_75_genome_counts'] == 1))]
    return c75_info


def get_pruned_c99(c99_info, opts='max', cutoff='0.4'):
    """ Prune centroids_99 shorter than 40% of the the max of corresponding clusters at 95 level """
    cluster_95_length = c99_info.groupby('centroid_95')['centroid_99_gene_length'].agg(lambda x: getattr(x, opts)()).reset_index()
    cluster_95_length.columns = ['centroid_95', 'cluster_length']
    c99_info = c99_info.merge(cluster_95_length, on='centroid_95')
    pruned_c99 = c99_info[c99_info['centroid_99_gene_length'] >= c99_info['cluster_length'] * float(cutoff)]
    return pruned_c99


def write_pruned_c99_seqs(from_fp, to_fp, c99s_to_keep):
    with open(to_fp, 'w') as ostream, \
        InputStream(from_fp, check_path=False) as istream:
        for rec in Bio.SeqIO.parse(istream, 'fasta'):
            gid = rec.id
            if gid in c99s_to_keep:
                gseq = str(rec.seq).upper()
                ostream.write(f">{gid}\n{gseq}\n")


def prune_centroids(args):
    if args.zzz_worker_mode:
        prune_centroids_worker(args)
    else:
        prune_centroids_master(args)


def prune_centroids_master(args):

    midas_db = MIDAS_DB(os.path.abspath(args.midasdb_dir), args.midasdb_name)
    species = midas_db.uhgg.species
    num_threads = args.num_threads

    def species_work(species_id):
        """
        For each species, we assume the user already ran the per genome annotate_pangenome.
        Therefore, here we just gather all the genomes' headerless file, and cat into one
        for each functional annotation method.
        """
        assert species_id in species, f"Species {species_id} is not in the database."
        species_genomes = species[species_id]
        msg = f"Pruning centorids_99 sequences for species {species_id} with {len(species_genomes)} total genomes."

        last_dest = midas_db.get_target_layout("pruned_centroids", False, species_id, args.prune_method, args.prune_cutoff)
        if os.path.exists(last_dest):
            if not args.force:
                tsprint(f"Destination {last_dest} for species {species_id} pangenome already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Pruning", "Repruning")
        tsprint(msg)

        local_dir = os.path.dirname(last_dest)
        command(f"mkdir -p {local_dir}")

        local_dest = midas_db.get_target_layout("prune_log", False, species_id)
        local_dir = os.path.dirname(local_dest)
        command(f"mkdir -p {local_dir}")

        worker_log = os.path.basename(local_dest)
        worker_subdir = os.path.dirname(local_dest) if args.scratch_dir == "." else f"{args.scratch_dir}"
        worker_log = f"{worker_subdir}/{worker_log}"

        if not args.debug:
            command(f"rm -rf {worker_subdir}")

        # Recurisve call via subcommand.  Use subdir, redirect logs.
        subcmd_str = f"--zzz_worker_mode -t {num_threads} --prune_method {args.prune_method} --prune_cutoff {args.prune_cutoff} {'--remove_singleton' if args.remove_singleton else ''} --midasdb_name {args.midasdb_name} --midasdb_dir {os.path.abspath(args.midasdb_dir)} {'--debug' if args.debug else ''} --scratch_dir {args.scratch_dir}"
        worker_cmd = f"cd {worker_subdir}; PYTHONPATH={pythonpath()} {sys.executable} -m midas2 prune_centroids -s {species_id} {subcmd_str} &>> {worker_log}"

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
    multithreading_map(species_work, species_id_list, args.num_threads)


def prune_centroids_worker(args):

    violation = "Please do not call prune_centroids_worker directly.  Violation"
    assert args.zzz_worker_mode, f"{violation}:  Missing --zzz_worker_mode arg."

    species_id = args.species
    midas_db = MIDAS_DB(args.midasdb_dir, args.midasdb_name)

    species = midas_db.uhgg.species
    assert species_id in species, f"{violation}: Species {species_id} is not in the database."

    c75_info = pd.read_csv(local_file(midas_db, species_id, "75"), sep="\t")
    if args.remove_singleton:
        c75_info = remove_singleton_c75_ids(c75_info)
    c75_ids_pass = set(c75_info['centroid_75'])

    c99_info = pd.read_csv(local_file(midas_db, species_id, "99"), sep="\t")
    c99_info = c99_info[['centroid_99', 'centroid_95', 'centroid_75', 'centroid_99_gene_length']]

    opts = args.prune_method
    cutoff = args.prune_cutoff
    pruned_c99 = get_pruned_c99(c99_info, opts, cutoff)
    pruned_c99 = pruned_c99[pruned_c99['centroid_75'].isin(c75_ids_pass)]
    c99_ids_keep = set(pruned_c99['centroid_99'])

    # Fetch centroids_99 sequences
    c99_ffn = midas_db.get_target_layout("pangenome_centroids", False, species_id)
    if_remove_singleton = 'rmsig.' if args.remove_singleton else ''
    local_dest = f"centroids_by.{opts}.{cutoff}.{if_remove_singleton}ffn"
    write_pruned_c99_seqs(c99_ffn, local_dest, c99_ids_keep)
    copy_tasks=[(local_dest, localpath(midas_db, species_id, f"pruned/{local_dest}"))]

    if args.scratch_dir != ".":
        multithreading_map(copy_star, copy_tasks, 4)


@register_args
def main(args):
    tsprint(f"Executing midas2 subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    prune_centroids(args)
