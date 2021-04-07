#!/usr/bin/env python3
import os
import sys
import Bio.SeqIO
from iggtools.common.argparser import add_subcommand, SUPPRESS
from iggtools.common.utils import tsprint, InputStream, retry, command, multithreading_map, find_files, upload, pythonpath, upload_star, num_physical_cores, download_reference
from iggtools.models.uhgg import UHGG, get_uhgg_layout, destpath
from iggtools.params import outputs
from iggtools.params.schemas import hmmsearch_max_evalue, hmmsearch_min_cov
from iggtools.subcommands.import_uhgg import decode_genomes_arg


CONCURRENT_MARKER_GENES_IDENTIFY = num_physical_cores


@retry
def find_files_with_retry(f):
    return find_files(f)


@retry
def fetch_genes(annotated_genes):
    """" Lookup of seq_id to sequence for PATRIC genes """
    gene_seqs = {}
    with InputStream(annotated_genes) as genes:
        for rec in Bio.SeqIO.parse(genes, 'fasta'):
            gene_seqs[rec.id] = str(rec.seq).upper()
    return gene_seqs


def hmmsearch(genome_id, species_id, marker_genes_hmm, num_threads=1):
    """ Performance HMM search using prokka annotated protein sequences """

    input_annotations = destpath(get_uhgg_layout(species_id, "faa", genome_id)["annotation_file"])
    annotated_genes = download_reference(input_annotations)

    hmmsearch_file = f"{genome_id}.hmmsearch"

    # Command
    if find_files(hmmsearch_file):
        # This only happens in debug mode, where we can use pre-existing file.
        tsprint(f"Found hmmsearch results for genome {genome_id} from prior run.")
    else:
        try:
            command(f"hmmsearch --noali --cpu {num_threads} --domtblout {hmmsearch_file} {marker_genes_hmm} {annotated_genes}")
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


def compute_marker_genes(genome_id, species_id, marker_genes_hmm):

    hmmsearch_file = hmmsearch(genome_id, species_id, marker_genes_hmm, num_threads=1)

    input_annotations = destpath(get_uhgg_layout(species_id, "ffn", genome_id)["annotation_file"])
    genes = fetch_genes(input_annotations)

    # Parse local hmmsearch file
    hmmsearch_seq = f"{genome_id}.markers.fa"
    hmmsearch_map = f"{genome_id}.markers.map"
    with open(hmmsearch_seq, "w") as o_seq, open(hmmsearch_map, "w") as o_map:
        for rec in find_hits(hmmsearch_file):
            marker_gene = genes[rec["query"]].upper()
            marker_info = [species_id, genome_id, rec["query"], len(marker_gene), rec["target"]]
            o_map.write('\t'.join(str(mi) for mi in marker_info) + '\n')
            o_seq.write('>%s\n%s\n' % (rec['query'], marker_gene))
    return {"hmmsearch": hmmsearch_file, "markers.fa": hmmsearch_seq, "markers.map": hmmsearch_map}


def infer_markers(args):
    if args.zzz_slave_toc:
        infer_markers_slave(args)
    else:
        infer_markers_master(args)


def infer_markers_master(args):

    # Fetch table of contents and marker genes HMM model from s3.
    # This will be read separately by each species build subcommand, so we make a local copy.
    local_toc = download_reference(outputs.genomes)
    marker_genes_hmm = download_reference(destpath(get_uhgg_layout("")["marker_genes_hmm"]))

    db = UHGG(local_toc)
    species_for_genome = db.genomes

    def genome_work(genome_id):
        assert genome_id in species_for_genome, f"Genome {genome_id} is not in the database."
        species_id = species_for_genome[genome_id]

        dest_file = destpath(get_uhgg_layout(species_id, "markers.map", genome_id)["marker_genes"])
        msg = f"Running HMMsearch for genome {genome_id} from species {species_id}."
        if find_files_with_retry(dest_file):
            if not args.force:
                tsprint(f"Destination {dest_file} for genome {genome_id} already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Running", "Rerunning")

        tsprint(msg)
        logfile = get_uhgg_layout(species_id, "", genome_id)["marker_genes_log"]
        slave_log = os.path.basename(logfile)
        slave_subdir = f"{species_id}__{genome_id}"
        if not args.debug:
            command(f"rm -rf {slave_subdir}")
        if not os.path.isdir(slave_subdir):
            command(f"mkdir {slave_subdir}")

        # Recurisve call via subcommand.  Use subdir, redirect logs.
        slave_cmd = f"cd {slave_subdir}; PYTHONPATH={pythonpath()} {sys.executable} -m iggtools infer_markers --genome {genome_id} --zzz_slave_mode --zzz_slave_toc {os.path.abspath(local_toc)} --zzz_slave_marker_genes_hmm {os.path.abspath(marker_genes_hmm)} {'--debug' if args.debug else ''} &>> {slave_log}"
        with open(f"{slave_subdir}/{slave_log}", "w") as slog:
            slog.write(msg + "\n")
            slog.write(slave_cmd + "\n")
        try:
            command(slave_cmd)
        finally:
            # Cleanup should not raise exceptions of its own, so as not to interfere with any
            # prior exceptions that may be more informative.  Hence check=False.
            upload(f"{slave_subdir}/{slave_log}", destpath(logfile), check=False)
            if not args.debug:
                command(f"rm -rf {slave_subdir}", check=False)

    genome_id_list = decode_genomes_arg(args, species_for_genome)
    multithreading_map(genome_work, genome_id_list, num_threads=CONCURRENT_MARKER_GENES_IDENTIFY)


def infer_markers_slave(args):
    """
    https://github.com/czbiohub/iggtools/wiki
    """

    violation = "Please do not call build_merker_genes_slave directly.  Violation"
    assert args.zzz_slave_mode, f"{violation}:  Missing --zzz_slave_mode arg."
    assert os.path.isfile(args.zzz_slave_toc), f"{violation}: File does not exist: {args.zzz_slave_toc}"
    assert os.path.isfile(args.zzz_slave_marker_genes_hmm), f"{violation}: Maker genes HMM model file does not exist: {args.zzz_slave_marker_genes_hmm}"

    db = UHGG(args.zzz_slave_toc)
    species_for_genome = db.genomes

    genome_id = args.genomes
    species_id = species_for_genome[genome_id]
    marker_genes_hmm = args.zzz_slave_marker_genes_hmm

    output_files = compute_marker_genes(genome_id, species_id, marker_genes_hmm)

    # Upload to S3
    last_dest = destpath(get_uhgg_layout(species_id, "markers.map", genome_id)["marker_genes"])
    command(f"aws s3 rm --recursive {os.path.dirname(last_dest)}")

    upload_tasks = []
    for k, o in output_files.items():
        if k == "markers.map":
            continue
        upload_tasks.append((o, destpath(get_uhgg_layout(species_id, k, genome_id)["marker_genes"])))
    multithreading_map(upload_star, upload_tasks)

    # Upload this last because it indicates all other work has succeeded.
    upload(output_files["markers.map"], last_dest)


def register_args(main_func):
    subparser = add_subcommand('infer_markers', main_func, help='HMM search marker genes for given genomes')
    subparser.add_argument('--genomes',
                           dest='genomes',
                           required=False,
                           help="genome[,genome...] to import;  alternatively, slice in format idx:modulus, e.g. 1:30, meaning import genomes whose ids are 1 mod 30; or, the special keyword 'all' meaning all genomes")
    subparser.add_argument('--zzz_slave_toc',
                           dest='zzz_slave_toc',
                           required=False,
                           help=SUPPRESS) # "reserved to pass table of contents from master to slave"
    subparser.add_argument('--zzz_slave_marker_genes_hmm',
                           dest='zzz_slave_marker_genes_hmm',
                           required=False,
                           help=SUPPRESS) # "reserved to common database from master to slave"
    return main_func


@register_args
def main(args):
    tsprint(f"Executing iggtools subcommand {args.subcommand} with args {vars(args)}.")
    infer_markers(args)
