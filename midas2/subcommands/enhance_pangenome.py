#!/usr/bin/env python3
import os
import sys
import pandas as pd

from midas2.common.argparser import add_subcommand
from midas2.common.utils import tsprint, command, OutputStream, multithreading_map, cat_files, pythonpath, num_physical_cores, copy_star
from midas2.common.utilities import decode_species_arg, decorate_genes_info_with_annot, annotation_ratio_x_members
from midas2.models.midasdb import MIDAS_DB
from midas2.params.inputs import MIDASDB_NAMES
from midas2.params.schemas import COLS_GENOMAD, COLS_MEFINDER, COLS_RESFINDER, COLS_EGGNOG


"""
For converting the raw annotation results into desirable TSV, there are two levels tasks.

We further augmented the genes_info.tsv with ratio of phage genes, plasmids genes, MGE genes, AMR genes, and eggNOG annotation.
- clusters_xx_info.tsv: deluxe info file
- genes_annotated.tsv: gene_is_phage, gene_is_plasmid, gene_is_amr, genes_is_me

For species (--species), we cat all the above produced temp files into one file, and add header.
    - annotation/{species_id}/genomad_virus.tsv
    - annotation/{species_id}/genomad_plasmid.tsv
    - annotation/{species_id}/mefinder.tsv
    - annotation/{species_id}/resfinder.tsv
    - annotation/{species_id}/eggnog.tsv
    - annotation/cluster_xx_annot:tsv: 'phage_ratio', 'plasmid_ratio', 'amr_ratio', 'me_ratio', basic eggNOG annotation columns
"""


def localpath(midas_db, species_id, mstep):
    return midas_db.get_target_layout(f"panannot_{mstep}", False, species_id)


def cat_annotations_x_genomes(midas_db, species_id, annotation_tools):
    # Get the list of all genomes of given species
    list_of_genomes = set(midas_db.uhgg.species[species_id].keys())
    for mstep in annotation_tools:
        list_of_temp_files = [midas_db.get_target_layout("panannot_tempfile", False, species_id, gid, mstep) for gid in list_of_genomes]
        assert all(os.path.exists(fp) for fp in list_of_temp_files), "Missing genome level parsed files. Need to annotate_pangenome."

        last_file = f"{mstep}.tsv"
        with OutputStream(last_file) as stream:
            if "genomad_virus" == mstep:
                cols = COLS_GENOMAD
            if "genomad_plasmid" == mstep:
                cols = COLS_GENOMAD
            if "mefinder" == mstep:
                cols = COLS_MEFINDER
            if "resfinder" == mstep:
                cols = COLS_RESFINDER
            if "eggnog" == mstep:
                cols = COLS_EGGNOG
            stream.write('\t'.join(cols) + '\n')
        cat_files(list_of_temp_files, last_file, 8)


def enhance_genes_info(midas_db, species_id):
    # cat temp/genome_annotation_files into annotation/*.tsv
    annotation_tools = ['genomad_virus', 'genomad_plasmid', 'mefinder', 'resfinder', 'eggnog']
    cat_annotations_x_genomes(midas_db, species_id, annotation_tools)

    copy_tasks = []
    for mstep in annotation_tools:
        copy_tasks.append((f"{mstep}.tsv", localpath(midas_db, species_id, mstep)))

    # Decorate genes_info with binary gene annotation
    genes_info_fp = midas_db.get_target_layout("pangenome_genes_info", False, species_id)
    genes_info = pd.read_csv(genes_info_fp, sep="\t")
    annot_files = {mstep:localpath(midas_db, species_id, mstep) for mstep in annotation_tools[:4]}

    genes_decorated = decorate_genes_info_with_annot(genes_info, annot_files)
    return genes_decorated, copy_tasks


def enhance_pangenome(args):
    if args.zzz_worker_mode:
        enhance_pangenome_worker(args)
    else:
        enhance_pangenome_master(args)


def enhance_pangenome_master(args):

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
        msg = f"Enhancing pangenome for species {species_id} with {len(species_genomes)} total genomes."

        dest_file = midas_db.get_target_layout("pangenome_genes_annot", True, species_id)
        local_file = midas_db.get_target_layout("pangenome_genes_annot", False, species_id)
        local_dir = os.path.dirname(local_file)
        command(f"mkdir -p {local_dir}/annotation")

        if os.path.exists(dest_file):
            if not args.force:
                tsprint(f"Destination {dest_file} for species {species_id} pangenome already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Enhancing", "Reenhancing")
        if not args.upload and os.path.exists(local_file):
            if not args.force:
                tsprint(f"Destination {local_file} for species {species_id} genes_annotated already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Enhancing", "Reenhancing")

        tsprint(msg)
        last_dest = midas_db.get_target_layout("enhance_log", True, species_id)
        local_dest = midas_db.get_target_layout("enhance_log", False, species_id)
        local_dir = os.path.dirname(local_dest)
        command(f"mkdir -p {local_dir}/annotation")

        worker_log = os.path.basename(local_dest)
        worker_subdir = os.path.dirname(local_dest) if args.scratch_dir == "." else f"{args.scratch_dir}"
        worker_log = f"{worker_subdir}/{worker_log}"
        if not args.debug:
            command(f"rm -rf {worker_subdir}")

        worker_subdir = f"{worker_subdir}/annotation"
        if not os.path.isdir(f"{worker_subdir}"):
            command(f"mkdir -p {worker_subdir}")

        # Recurisve call via subcommand.  Use subdir, redirect logs.
        subcmd_str = f"--zzz_worker_mode -t {num_threads} --midasdb_name {args.midasdb_name} --midasdb_dir {os.path.abspath(args.midasdb_dir)} {'--debug' if args.debug else ''} {'--upload' if args.upload else ''} --scratch_dir {args.scratch_dir}"
        worker_cmd = f"cd {worker_subdir}; PYTHONPATH={pythonpath()} {sys.executable} -m midas2 enhance_pangenome -s {species_id} {subcmd_str} &>> {worker_log}"

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


def enhance_pangenome_worker(args):

    violation = "Please do not call enhance_pangenome_worker directly.  Violation"
    assert args.zzz_worker_mode, f"{violation}:  Missing --zzz_worker_mode arg."

    species_id = args.species
    midas_db = MIDAS_DB(args.midasdb_dir, args.midasdb_name, num_cores=1)

    species = midas_db.uhgg.species
    assert species_id in species, f"{violation}: Species {species_id} is not in the database."

    genes_annotated, copy_tasks = enhance_genes_info(midas_db, species_id)
    genes_annotated.to_csv("genes_annotated.tsv", sep='\t', index=False)

    if args.scratch_dir == ".":
        copy_tasks = []
    copy_tasks.append(("genes_annotated.tsv", midas_db.get_target_layout("pangenome_genes_annot", False, species_id)))

    # Compute centroid_xx annotationdensity
    for xx in ["99", "95", "90", "85", "80", "75"]:
        cxx_df = annotation_ratio_x_members(genes_annotated, localpath(midas_db, species_id, "eggnog"), xx)
        cxx_df.to_csv(f"clusters_{xx}_annot.tsv",  sep='\t', index=False)
        if args.scratch_dir != ".":
            copy_tasks.append((f"clusters_{xx}_annot.tsv", midas_db.get_target_layout("cluster_xx_annot", False, species_id, "", xx)))

    # Merge the cluster info with annotation
    for xx in ["99", "95", "90", "85", "80", "75"]:
        xx_info = pd.read_csv(midas_db.get_target_layout("cluster_xx_info", False, species_id, "", xx), sep="\t")
        xx_anno = pd.read_csv(midas_db.get_target_layout("cluster_xx_annot", False, species_id, "", xx), sep="\t")
        assert xx_info.shape[0] == xx_anno.shape[0], f"The shape of augment/clusters_{xx}_info.tsv disagree with annotation/clusters_{xx}_annoted.tsv"
        xx_df = pd.merge(xx_info, xx_anno, on=f"centroid_{xx}", how='left')
        # Write directly to target output
        xx_df.to_csv(midas_db.get_target_layout("pangenome_cluster_xx", False, species_id, "", xx), sep='\t', index=False)

    multithreading_map(copy_star, copy_tasks, 4)
    command("rm -f genes_annotated.tsv", check=False)


def register_args(main_func):
    subparser = add_subcommand('enhance_pangenome', main_func, help='Enhance pangenome for specified species')
    subparser.add_argument('-s',
                           '--species',
                           dest='species',
                           required=False,
                           help="species[,species...] whose pangenome(s) to build;  alternatively, species slice in format idx:modulus, e.g. 1:30, meaning enhance species whose ids are 1 mod 30; or, the special keyword 'all' meaning all species")
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
    tsprint(f"Executing midas2 subcommand {args.subcommand}.")
    enhance_pangenome(args)
