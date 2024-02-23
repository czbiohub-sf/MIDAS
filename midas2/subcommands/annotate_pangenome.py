#!/usr/bin/env python3
import os
import sys
from itertools import chain
import pandas as pd
from pybedtools import BedTool

from midas2.common.argparser import add_subcommand
from midas2.common.utils import tsprint, command, multithreading_map, pythonpath, num_physical_cores, copy_star
from midas2.common.utilities import decode_genomes_arg, decode_species_arg, scan_eggnog
from midas2.models.midasdb import MIDAS_DB
from midas2.params.inputs import MIDASDB_NAMES


"""
For converting the raw annotation results into desirable TSV, there are two levels tasks.
For genomes (--genomes), we sequentially parse the raw outputs of four MGE tools, and write parsed results to temp/.
"""


def make_subdirs(local_dir):
    command(f"mkdir -p {local_dir}/resfinder")
    command(f"mkdir -p {local_dir}/mefinder")
    command(f"mkdir -p {local_dir}/genomad_virus")
    command(f"mkdir -p {local_dir}/genomad_plasmid")
    command(f"mkdir -p {local_dir}/eggnog")


def parse_genomad_virus_genes(file_path):
    """ Parse Genomad virus genes TSV into DataFrame """
    df = pd.read_csv(file_path, delimiter="\t")

    if df.empty:
        return pd.DataFrame()

    df = df[['gene', 'start', 'end', 'annotation_conjscan', 'annotation_amr', 'annotation_accessions', 'annotation_description']]
    df[['contig_plus', 'gene_num']] = df['gene'].str.rsplit('_', n=1, expand=True)

    # For the provisus case, we only keep the contig_id
    if df['contig_plus'].str.contains(r'\|p').any():
        split_data = df['contig_plus'].str.split(r'\|p', expand=True)
        # If the delimiter is not present, then the second column will be NaN.
        # We'll use that to replace NaNs in the first column with the original 'contig_plus' values.
        split_data[0].where(split_data[1].notna(), df['contig_plus'], inplace=True)
        # Assign the results back to the original dataframe
        df['contig_external'] = split_data[0]
    else:
        df['contig_external'] = df['contig_plus']

    assert ~df[['contig_external', 'start', 'end']].duplicated().any(), f"Duplicated virus results for {file_path}"

    df = df[['contig_external', 'start', 'end', 'gene',  'annotation_conjscan', 'annotation_amr', 'annotation_accessions', 'annotation_description']]
    df = df.rename(columns={'start': 'start_genomad', 'end': 'end_genomad', 'gene':'gene_genomad'})

    return df


def parse_genomad_plasmid_genes(file_path):
    """ Parse Genomad plasmid genes TSV into DataFrame """
    df = pd.read_csv(file_path, delimiter="\t")

    if df.empty:
        return pd.DataFrame()

    df = df[['gene', 'start', 'end', 'annotation_conjscan', 'annotation_amr', 'annotation_accessions', 'annotation_description']]
    df[['contig_external', 'gene_num']] = df['gene'].str.rsplit('_', n=1, expand=True)

    assert ~df[['contig_external', 'start', 'end']].duplicated().any(), f"Duplicated plasmid results for {file_path}"

    df = df[['contig_external', 'start', 'end', 'gene', 'annotation_conjscan', 'annotation_amr', 'annotation_accessions', 'annotation_description']]
    df = df.rename(columns={'start': 'start_genomad', 'end': 'end_genomad', 'gene': 'gene_genomad'})

    return df


def parse_mefinder(file_path):
    """ Parse MEFINDER into DataFrame """
    df = pd.read_csv(file_path, delimiter=",", skiprows=5)

    if df.empty:
        return pd.DataFrame()

    df['contig_id'] = df['contig'].str.split(' ').str[0]
    df = df[['contig_id', 'start', 'end', 'mge_no', 'prediction', 'name', 'type', 'synonyms']]
    df = df.rename(columns={'contig_id': 'contig_external', 'start': 'start_mefinder', 'end': 'end_mefinder'})

    return df


def parse_resfinder(file_path):
    # rename the columns without spaces
    new_columns = ['resistance_gene', 'identity', 'align_ratio', 'coverage', 'within_reference_position', 'contig', 'within_contig_position', 'phenotype', 'accession_no']
    df = pd.read_csv(file_path, sep='\t', header=0, names=new_columns)

    if df.empty:
        return pd.DataFrame()

    df['contig'] = df['contig'].str.split(' ').str[0]
    df[['start', 'end']]  = df['within_contig_position'].str.split('\\.\\.', expand=True)
    # Notes there can be duplicated annotaitons for the same gene range (contig-start-end)
    df = df[['contig', 'start', 'end', 'resistance_gene', 'phenotype', 'accession_no']]
    df = df.rename(columns={'contig': 'contig_external', 'start': 'start_resfinder', 'end': 'end_resfinder'})
    return df


def merge_annot_with_genes(df, all_genes):
    bed2 = BedTool.from_dataframe(df)
    bed1 = BedTool.from_dataframe(all_genes)
    overlaps = bed1.intersect(bed2, wa=True, wb=True)

    if len(overlaps) == 0:
        return pd.DataFrame()

    overlapping_df = pd.read_table(overlaps.fn, header=None)
    overlapping_df.columns = list(all_genes.columns) + list(df.columns)
    # Drop duplicated column: contig_id
    overlapping_df = overlapping_df.drop(columns=['contig_external'])
    # Reorder the column names
    overlapping_df = overlapping_df[['gene_id'] + [col for col in overlapping_df if col != 'gene_id']]
    return overlapping_df


def write_processed_genome_annot(df, all_genes, local_dest):
    if df.empty:
        command(f"touch {local_dest}")
    else:
        df = merge_annot_with_genes(df, all_genes)
        if df.empty:
            command(f"touch {local_dest}")
        else:
            df.to_csv(local_dest, sep='\t', index=False, header=False)


def annotate_pangenome(args):
    if args.zzz_worker_mode:
        annotate_pangenome_worker(args)
    else:
        annotate_pangenome_master(args)


def annotate_pangenome_master(args):

    midas_db = MIDAS_DB(os.path.abspath(args.midasdb_dir), args.midasdb_name)
    species = midas_db.uhgg.species
    species_for_genome = midas_db.uhgg.genomes

    def genome_work(genome_id):
        """
        For each genome, we sequentially read in three functional annotation tables,
        overlap with the list_of_centroids99, write to TEMP files.
        """

        assert genome_id in species_for_genome, f"Genome {genome_id} is not in the database."
        species_id = species_for_genome[genome_id]

        local_dest = midas_db.get_target_layout("panannot_tempfile", False, species_id, genome_id, "genomad_virus")
        local_dir = os.path.dirname(os.path.dirname(local_dest))

        msg = f"Parsing genome {genome_id} from species {species_id}."
        if os.path.exists(local_dest):
            if not args.force:
                tsprint(f"Destination {local_dest} for genome {genome_id} annotations already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Parsing", "Reparsing")
        tsprint(msg)

        if not args.debug:
            command(f"rm -rf {local_dir}")
        make_subdirs(local_dir)

        worker_subdir = local_dir if args.scratch_dir == "." else f"{args.scratch_dir}/functional_annotation/{species_id}"
        if not args.debug:
            command(f"rm -rf {worker_subdir}")
        if not os.path.isdir(worker_subdir):
            command(f"mkdir -p {worker_subdir}")
        make_subdirs(worker_subdir)

        # Recurisve call via subcommand.  Use subdir, redirect logs.
        # Output files are generated inside worker_subdir
        worker_log = f"{worker_subdir}/{genome_id}.log"
        subcmd_str = f"--zzz_worker_mode --midasdb_name {args.midasdb_name} --midasdb_dir {os.path.abspath(args.midasdb_dir)} {'--debug' if args.debug else ''} --scratch_dir {args.scratch_dir}"
        worker_cmd = f"cd {worker_subdir}; PYTHONPATH={pythonpath()} {sys.executable} -m midas2 annotate_pangenome --genome {genome_id} {subcmd_str} &>> {worker_log}"

        try:
            command(worker_cmd)
        finally:
            command(f"rm -f {worker_log}")

    if args.species:
        species_id_list = decode_species_arg(args, species)
        genome_id_list = list(chain.from_iterable(species[spid].keys() for spid in species_id_list))
    if args.genomes:
        genome_id_list = decode_genomes_arg(args, species_for_genome)
    multithreading_map(genome_work, genome_id_list, num_threads=args.num_threads)


def annotate_pangenome_worker(args):
    """
    For given genome, we take the following inputs:
    Input:
        - Species: pangenomes/cluster_info.txt, contigs_len
        - Genome: genes.feature
        - Funannot: genomad, mefinder and resfinder results
    Output:
        - Four headerless temp tables for centroids.99
    """

    violation = "Please do not call annotate_pangenome_worker directly.  Violation"
    assert args.zzz_worker_mode, f"{violation}:  Missing --zzz_worker_mode arg."

    midas_db = MIDAS_DB(args.midasdb_dir, args.midasdb_name)
    species_for_genome = midas_db.uhgg.genomes

    genome_id = args.genomes
    species_id = species_for_genome[genome_id]

    cluster_info_fp = midas_db.get_target_layout("cluster_xx_info", False, species_id, "", "99")
    contig_len_fp = midas_db.get_target_layout("pangenome_contigs_len", False, species_id)
    gene_feature_fp = midas_db.get_target_layout("annotation_genes", False, species_id, genome_id)

    contig_len = pd.read_csv(contig_len_fp, sep='\t')
    gene_features = pd.read_csv(gene_feature_fp, sep='\t')

    # 2023-12-05: we keep ALL the GENES, instead of centroids_99
    all_genes = pd.merge(gene_features, contig_len[['contig_id', 'contig_length']], left_on='contig_id', right_on='contig_id', how='inner')
    all_genes = all_genes[['contig_id', 'start', 'end', 'gene_id', 'strand', 'gene_type', 'contig_length']]

    # Genomad virus
    df = parse_genomad_virus_genes(midas_db.get_target_layout("genomad_virus_genes", False, species_id, genome_id))
    local_dest = f"genomad_virus/{genome_id}"
    write_processed_genome_annot(df, all_genes, local_dest)

    # Genomad plasmid
    df = parse_genomad_plasmid_genes(midas_db.get_target_layout("genomad_plasmid_genes", False, species_id, genome_id))
    local_dest = f"genomad_plasmid/{genome_id}"
    write_processed_genome_annot(df, all_genes, local_dest)

    # MEfinder
    df = parse_mefinder(midas_db.get_target_layout("mefinder_results", False, species_id, genome_id))
    local_dest = f"mefinder/{genome_id}"
    write_processed_genome_annot(df, all_genes, local_dest)

    # Resfinder
    df = parse_resfinder(midas_db.get_target_layout("resfinder_results", False, species_id, genome_id))
    local_dest = f"resfinder/{genome_id}"
    write_processed_genome_annot(df, all_genes, local_dest)

    # EggNOG was ran on the centroid_99s.
    centroids_99 = pd.read_csv(cluster_info_fp, sep='\t', usecols=[0]) # Only keep the list of centroid_99
    centroids_99 = pd.merge(centroids_99, gene_features, left_on='centroid_99', right_on='gene_id', how='inner')
    centroids_99 = pd.merge(centroids_99, contig_len[['contig_id', 'contig_length']], left_on='contig_id', right_on='contig_id', how='inner')
    centroids_99 = centroids_99[['contig_id', 'start', 'end', 'centroid_99', 'strand', 'gene_type', 'contig_length']]

    df = scan_eggnog(midas_db.get_target_layout("eggnog_results", False, species_id))
    local_dest = f"eggnog/{genome_id}"
    merged_df  = df.merge(centroids_99, left_on='#query', right_on='centroid_99', how='inner')
    merged_df = merged_df.drop(columns=['centroid_99'])
    merged_df.to_csv(local_dest, sep='\t', index=False, header=False)

    if args.scratch_dir != ".":
        copy_tasks = [
            (f"genomad_virus/{genome_id}", midas_db.get_target_layout("panannot_tempfile", False, species_id, genome_id, "genomad_virus")),
            (f"genomad_plasmid/{genome_id}", midas_db.get_target_layout("panannot_tempfile", False, species_id, genome_id, "genomad_plasmid")),
            (f"mefinder/{genome_id}", midas_db.get_target_layout("panannot_tempfile", False, species_id, genome_id, "mefinder")),
            (f"resfinder/{genome_id}", midas_db.get_target_layout("panannot_tempfile", False, species_id, genome_id, "resfinder")),
            (f"eggnog/{genome_id}", midas_db.get_target_layout("panannot_tempfile", False, species_id, genome_id, "eggnog"))
        ]
        multithreading_map(copy_star, copy_tasks, 2)


def register_args(main_func):
    subparser = add_subcommand('annotate_pangenome', main_func, help='Decorate pangenome with annotated genomes.')
    subparser.add_argument('-s',
                           '--species',
                           dest='species',
                           required=False,
                           help="species[,species...] whose pangenome(s) to build;  alternatively, species slice in format idx:modulus, e.g. 1:30, meaning build species whose ids are 1 mod 30; or, the special keyword 'all' meaning all species")
    subparser.add_argument('--genomes',
                           dest='genomes',
                           required=False,
                           help="genome[,genome...] to import;  alternatively, slice in format idx:modulus, e.g. 1:30, meaning pangenome annotate whose ids are 1 mod 30; or, the special keyword 'all' meaning all genomes")
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
    return main_func


@register_args
def main(args):
    tsprint(f"Executing midas2 subcommand {args.subcommand}.")
    annotate_pangenome(args)
