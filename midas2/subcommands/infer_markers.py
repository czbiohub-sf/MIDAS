#!/usr/bin/env python3
import os
import sys
from midas2.common.argparser import add_subcommand, SUPPRESS
from midas2.common.utils import tsprint, InputStream, retry, command, multithreading_map, find_files, upload, pythonpath, upload_star, num_physical_cores
from midas2.common.utilities import scan_genes, decode_genomes_arg
from midas2.models.midasdb import MIDAS_DB
from midas2.params.inputs import hmmsearch_max_evalue, hmmsearch_min_cov, MIDASDB_NAMES


CONCURRENT_INFER_MARKERS = num_physical_cores


@retry
def find_files_with_retry(f):
    return find_files(f)


def hmm_search(genome_id, input_annotations, marker_genes_hmm, num_threads=1):
    """ Performance HMM search using prokka annotated protein sequences """
    hmmsearch_file = f"{genome_id}.hmmsearch"
    if find_files(hmmsearch_file):
        # This only happens in debug mode, where we can use pre-existing file.
        tsprint(f"Found hmmsearch results for genome {genome_id} from prior run.")
    else:
        try:
            command(f"hmmsearch --noali --cpu {num_threads} --domtblout {hmmsearch_file} {marker_genes_hmm} {input_annotations}")
        except:
            # Do not keep bogus zero-length files;  those are harmful if we rerun in place.
            command(f"mv {hmmsearch_file} {hmmsearch_file}.bogus", check=False)
            raise
    return hmmsearch_file


def parse_hmmsearch(hmmsearch_file):
    """ Parse HMMER domblout files. Return data-type formatted dictionary """
    with InputStream(hmmsearch_file) as f_in:
        for line in f_in:
            if line[0] == "#":
                continue
            x = line.rstrip().split()
            query = x[0]
            target = x[3]
            evalue = float(x[12])
            qcov = (int(x[20]) - int(x[19]) + 1)/float(x[2])
            tcov = (int(x[16]) - int(x[15]) + 1)/float(x[5])
            yield {'query':query, 'target':target, 'evalue':evalue, 'qcov':qcov, 'tcov':tcov, 'qlen':int(x[2]), 'tlen':int(x[5])}


def find_hits(hmmsearch_file):
    hits = {}
    for r in parse_hmmsearch(hmmsearch_file):
        if r['evalue'] > hmmsearch_max_evalue:
            continue
        if min(r['qcov'], r['tcov']) < hmmsearch_min_cov:
            continue
        if r['target'] not in hits:
            hits[r['target']] = r
        elif r['evalue'] < hits[r['target']]['evalue']:
            hits[r['target']] = r
    return list(hits.values())


def compute_marker_genes(genome_id, species_id, marker_genes_hmm, midas_db):
    """ Parse local HMM search output file """
    s3_gene_faa = midas_db.fetch_file("annotation_faa", species_id, genome_id)
    hmmsearch_file = hmm_search(genome_id, s3_gene_faa, marker_genes_hmm, num_threads=1)

    s3_gene_seq = midas_db.fetch_file("annotation_ffn", species_id, genome_id)
    genes = scan_genes(s3_gene_seq)

    hmmsearch_seq = f"{genome_id}.markers.fa"
    hmmsearch_map = f"{genome_id}.markers.map"
    with open(hmmsearch_seq, "w") as o_seq, open(hmmsearch_map, "w") as o_map:
        for rec in find_hits(hmmsearch_file):
            marker_gene = genes[rec["query"]].upper()
            marker_info = [species_id, genome_id, rec["query"], len(marker_gene), rec["target"]]
            o_map.write('\t'.join(str(mi) for mi in marker_info) + '\n')
            o_seq.write('>%s\n%s\n' % (rec['query'], marker_gene))
    return {"marker_genes_hmmsearch": hmmsearch_file, "marker_genes_seq": hmmsearch_seq, "marker_genes_map": hmmsearch_map}


def infer_markers(args):
    if args.zzz_worker_mode:
        infer_markers_worker(args)
    else:
        infer_markers_master(args)


def infer_markers_master(args):

    # Fetch table of contents and marker genes HMM model from s3.
    # This will be read separately by each species build subcommand, so we make a local copy.
    midas_db = MIDAS_DB(os.path.abspath(args.midasdb_dir), args.midasdb_name)
    species_for_genome = midas_db.uhgg.genomes

    marker_genes_hmm = midas_db.fetch_files("marker_db_hmm")

    def genome_work(genome_id):
        assert genome_id in species_for_genome, f"Genome {genome_id} is not in the database."
        species_id = species_for_genome[genome_id]

        dest_file = midas_db.get_target_layout("marker_genes_map", True, species_id, genome_id)

        msg = f"Running HMMsearch for genome {genome_id} from species {species_id}."
        if find_files_with_retry(dest_file):
            if not args.force:
                tsprint(f"Destination {dest_file} for genome {genome_id} already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Running", "Rerunning")

        tsprint(msg)
        last_dest = midas_db.get_target_layout("marker_genes_log", True, species_id, genome_id)

        worker_log = midas_db.get_target_layout("marker_genes_log", False, species_id, genome_id)
        worker_subdir = os.path.dirname(worker_log)

        if not args.debug:
            command(f"rm -rf {worker_subdir}")
        if not os.path.isdir(worker_subdir):
            command(f"mkdir -p {worker_subdir}")

        # Recurisve call via subcommand.  Use subdir, redirect logs.
        worker_cmd = f"cd {worker_subdir}; PYTHONPATH={pythonpath()} {sys.executable} -m midas2 infer_markers --genome {genome_id} --midasdb_name {args.midasdb_name} --midasdb_dir {os.path.abspath(args.midasdb_dir)} --zzz_worker_mode --zzz_worker_marker_genes_hmm {os.path.abspath(marker_genes_hmm)} {'--debug' if args.debug else ''} &>> {worker_log}"
        with open(f"{worker_log}", "w") as slog:
            slog.write(msg + "\n")
            slog.write(worker_cmd + "\n")

        try:
            command(worker_cmd)
        finally:
            # Cleanup should not raise exceptions of its own, so as not to interfere with any
            # prior exceptions that may be more informative.  Hence check=False.
            if not args.debug:
                upload(f"{worker_log}", last_dest, check=False)
                command(f"rm -rf {worker_subdir}", check=False)

    genome_id_list = decode_genomes_arg(args, species_for_genome)
    multithreading_map(genome_work, genome_id_list, num_threads=CONCURRENT_INFER_MARKERS)


def infer_markers_worker(args):
    """
    https://github.com/czbiohub/MIDAS2.0/wiki/MIDAS-DB
    """

    violation = "Please do not call infer_markers_worker directly.  Violation"
    assert args.zzz_worker_mode, f"{violation}:  Missing --zzz_worker_mode arg."
    assert os.path.isfile(args.zzz_worker_marker_genes_hmm), f"{violation}: Maker genes HMM model file does not exist: {args.zzz_worker_marker_genes_hmm}"

    midas_db = MIDAS_DB(args.midasdb_dir, args.midasdb_name)
    species_for_genome = midas_db.uhgg.genomes

    genome_id = args.genomes
    species_id = species_for_genome[genome_id]
    marker_genes_hmm = args.zzz_worker_marker_genes_hmm

    output_files = compute_marker_genes(genome_id, species_id, marker_genes_hmm, midas_db)

    # Upload to S3
    if not args.debug:
        last_dest = midas_db.get_target_layout("marker_genes_map", True, species_id, genome_id)
        command(f"aws s3 rm --recursive {os.path.dirname(last_dest)}")

        upload_tasks = []
        for k, o in output_files.items():
            if k == "marker_genes_map":
                continue
            upload_tasks.append((o, midas_db.get_target_layout(k, True, species_id, genome_id)))
        multithreading_map(upload_star, upload_tasks)

        # Upload this last because it indicates all other work has succeeded.
        upload(output_files["marker_genes_map"], last_dest)


def register_args(main_func):
    subparser = add_subcommand('infer_markers', main_func, help='HMM search marker genes for specified genomes')
    subparser.add_argument('--genomes',
                           dest='genomes',
                           required=False,
                           help="genome[,genome...] to import;  alternatively, slice in format idx:modulus, e.g. 1:30, meaning import genomes whose ids are 1 mod 30; or, the special keyword 'all' meaning all genomes")
    subparser.add_argument('--zzz_worker_marker_genes_hmm',
                           dest='zzz_worker_marker_genes_hmm',
                           required=False,
                           help=SUPPRESS) # "reserved to common database from master to worker"
    subparser.add_argument('--midasdb_name',
                           dest='midasdb_name',
                           type=str,
                           default="uhgg",
                           choices=MIDASDB_NAMES,
                           help=f"MIDAS Database name.")
    subparser.add_argument('--midasdb_dir',
                           dest='midasdb_dir',
                           type=str,
                           default=".",
                           help=f"Local MIDAS Database path mirroing S3.")
    return main_func


@register_args
def main(args):
    tsprint(f"Executing midas2 subcommand {args.subcommand} with args {vars(args)}.")
    infer_markers(args)
