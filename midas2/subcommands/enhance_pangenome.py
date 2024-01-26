#!/usr/bin/env python3
import os
import sys
from multiprocessing import Semaphore
import Bio.SeqIO
import pandas as pd

from midas2.common.argparser import add_subcommand
from midas2.common.utils import InputStream, tsprint, command, OutputStream, multithreading_map, cat_files, pythonpath, num_physical_cores, copy_star
from midas2.common.utilities import decode_species_arg, decorate_genes_info, annotation_ratio_x_members
from midas2.models.midasdb import MIDAS_DB
from midas2.params.inputs import MIDASDB_NAMES
from midas2.params.schemas import COLS_GENOMAD, COLS_MEFINDER, COLS_RESFINDER, COLS_EGGNOG


"""
For converting the raw annotation results into desirable TSV, there are two levels tasks.
For species (--species), we cat all the above produced temp files into one file, and add header.
    - pangenomes/annotation/{species_id}/genomad_virus.tsv
    - pangenomes/annotation/{species_id}/genomad_plasmid.tsv
    - pangenomes/annotation/{species_id}/mefinder.tsv
    - pangenomes/annotation/{species_id}/resfinder.tsv
    - pangenomes/annotation/{species_id}/eggnog.tsv

We further augmented the genes_info.tsv with ratio of phage genes, plasmids genes, MGE genes, AMR genes, and eggNOG annotation.
- pangenomes/genes_annotated.tsv: gene_is_phage, gene_is_plasmid, gene_is_amr, genes_is_me
- pangenomes/centroid_xx_annotated:tsv: 'phage_ratio', 'plasmid_ratio', 'amr_ratio', 'me_ratio', basic eggNOG annotation columns

We then generated the pruned/centroids.99.ffn
"""


CONCURRENT_SPECIES_BUILDS = Semaphore(8)

def localpath(midas_db, species_id, mstep):
    return midas_db.get_target_layout(f"panannot_{mstep}", False, species_id)


def cat_annotations_x_genomes(midas_db, list_of_genomes, species_id):
    for mstep in ['genomad_virus', 'genomad_plasmid', 'mefinder', 'resfinder', 'eggnog']:
        list_of_temp_files = [midas_db.get_target_layout("panannot_tempfile", False, species_id, gid, mstep) for gid in list_of_genomes]
        assert all(os.path.exists(fp) for fp in list_of_temp_files), "Missing genome level parsed files. Need to rerun with --genomes all option."

        last_file = f"annotation/{mstep}.tsv"
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
        cat_files(list_of_temp_files, last_file, 10)


def is_float(s):
    try:
        float(s)  # Try to convert the string to a float
        return True
    except ValueError:
        return False


def prune_c99s(cluster_info_fp, opts = "0.4"):
    """ Prune centroids_99 shorter than 40% of the corresponding centroids_95"""
    # What do we need here: clusters_info.tsv: c99 - c95 - c99_length
    cluster_info = pd.read_csv(cluster_info_fp, sep='\t') # Only keep the list of centroid_99
    if is_float(opts):
        # Filter out centroids_99 whose gene length is shorter than 40% of the corresponding centroids_95
        group_c95_length = cluster_info[['centroid_95', 'centroid_99_length']].drop_duplicates()
        group_c95_length['centroid_99_length'] = group_c95_length['centroid_99_length'] * float(opts)
    else:
        if opts not in ('median', 'max', 'mean'):
            raise ValueError("opts must be 'median', 'max', or 'mean'")
        # Filter out centroids_99 whose gene length is shorter than median of all centroids_99 members
        group_c95_length = cluster_info.groupby('centroid_95')['centroid_99_length'].agg(lambda x: getattr(x, opts)())
    group_c95_length.columns = ['centroid_95', 'grouped_length']

    cluster_info = cluster_info.merge(group_c95_length, on='centroid_95')
    pruned_cluster = cluster_info[cluster_info['centroid_99_length'] > cluster_info['grouped_length']]

    c99s_to_keep = set(pruned_cluster['centroid_99'])
    return c99s_to_keep


def write_pruned_c99s(orig_ffn_fp, pruned_ffn_fp, c99s_to_keep):
    with open(pruned_ffn_fp, 'w') as ostream, \
        InputStream(orig_ffn_fp, check_path=False) as istream:
        for rec in Bio.SeqIO.parse(istream, 'fasta'):
            gid = rec.id
            if gid in c99s_to_keep:
                gseq = str(rec.seq).upper()
                ostream.write(f">{gid}\n{gseq}\n")


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

        dest_file = midas_db.get_target_layout("pangenome_genes_annot", False, species_id)
        msg = f"Enhancing pangenome for species {species_id} with {len(species_genomes)} total genomes."

        if os.path.exists(dest_file):
            if not args.force:
                tsprint(f"Destination {dest_file} for species {species_id} pangenome already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Enhancing", "Reenhancing")

        with CONCURRENT_SPECIES_BUILDS:
            tsprint(msg)
            last_dest = midas_db.get_target_layout("pangenome_enhance_log", True, species_id)
            local_dest = midas_db.get_target_layout("pangenome_enhance_log", False, species_id)
            local_dir = os.path.dirname(local_dest)
            command(f"mkdir -p {local_dir}/annotation")
            command(f"mkdir -p {local_dir}/pruned")

            worker_log = os.path.basename(local_dest)
            worker_subdir = os.path.dirname(local_dest) if args.scratch_dir == "." else f"{args.scratch_dir}"
            worker_log = f"{worker_subdir}/{worker_log}"

            if not args.debug:
                command(f"rm -rf {worker_subdir}")
            if not os.path.isdir(f"{worker_subdir}/pruned"):
                command(f"mkdir -p {worker_subdir}/annotation")
                command(f"mkdir -p {worker_subdir}/pruned")

            # Recurisve call via subcommand.  Use subdir, redirect logs.
            subcmd_str = f"--zzz_worker_mode -t {num_threads} --prune_cutoff {args.prune_cutoff}  --midasdb_name {args.midasdb_name} --midasdb_dir {os.path.abspath(args.midasdb_dir)} {'--debug' if args.debug else ''} {'--upload' if args.upload else ''} --scratch_dir {args.scratch_dir}"
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
    multithreading_map(species_work, species_id_list, args.num_threads)


def enhance_pangenome_worker(args):

    violation = "Please do not call enhance_pangenome_worker directly.  Violation"
    assert args.zzz_worker_mode, f"{violation}:  Missing --zzz_worker_mode arg."

    species_id = args.species
    midas_db = MIDAS_DB(args.midasdb_dir, args.midasdb_name)

    species = midas_db.uhgg.species
    assert species_id in species, f"{violation}: Species {species_id} is not in the database."

    # cat temp/genome_annotation_files into one
    list_of_genomes = set(species[species_id].keys())
    cat_annotations_x_genomes(midas_db, list_of_genomes, species_id)

    copy_tasks = []
    for mstep in ['genomad_virus', 'genomad_plasmid', 'mefinder', 'resfinder', 'eggnog']:
        copy_tasks.append((f"annotation/{mstep}.tsv", localpath(midas_db, species_id, mstep)))

    # decprate genes_info with binary gene annotation
    genes_info_fp = midas_db.get_target_layout("pangenome_genes_info", False, species_id)
    genes_decorated = decorate_genes_info(midas_db, genes_info_fp, species_id)
    genes_decorated.to_csv("genes_annotated.tsv", sep='\t', index=False)

    genes_anno_fp = midas_db.get_target_layout("pangenome_genes_annot", False, species_id)
    copy_tasks.append(("genes_annotated.tsv", genes_anno_fp))

    # compute centroid_xx_annotation_density
    for xx in ['95', '80', '75']:
        cxx_df = annotation_ratio_x_members(midas_db, species_id, xx)
        cxx_df.to_csv(f"centroids.{xx}_annotated.tsv",  sep='\t', index=False)

        cluster_annoted_fp = midas_db.get_target_layout("pangenome_cluster_annot", False, species_id, '', xx)
        copy_tasks.append((f"centroids.{xx}_annotated.tsv", cluster_annoted_fp))

    # prune centroids.99
    opts = args.prune_cutoff
    cluster_info_fp = midas_db.get_target_layout("pangenome_cluster_info", False, species_id)
    c99s_to_keep = prune_c99s(cluster_info_fp, opts)

    centroids99_fp = midas_db.get_target_layout("pangenome_centroids", False, species_id)
    write_pruned_c99s(centroids99_fp, f"pruned/centroids_opts.{opts}.ffn", c99s_to_keep)
    last_file = midas_db.get_target_layout("pruned_centroids", False, species_id, "", opts)
    copy_tasks.append((f"pruned/centroids_opts.{opts}.ffn", last_file))

    if args.scratch_dir != ".":
        multithreading_map(copy_star, copy_tasks, 4)


def register_args(main_func):
    subparser = add_subcommand('enhance_pangenome', main_func, help='Genome annotation for specified genomes using Prokka with all cores')
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
    subparser.add_argument('--prune_cutoff',
                           dest='prune_cutoff',
                           type=float,
                           default=0.4,
                           help="Prune cutoff: centroid99 shorter than 0.4 of the centroid95 are pruned for reading mapping.")
    return main_func


@register_args
def main(args):
    tsprint(f"Executing midas2 subcommand {args.subcommand}.")
    enhance_pangenome(args)
