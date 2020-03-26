import json
import os
import random
from collections import defaultdict
from itertools import chain
import numpy as np

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, num_physical_cores, InputStream, OutputStream, select_from_tsv, download_reference, command
from iggtools.models.uhgg import UHGG, fetch_marker_genes
from iggtools.models.sample import Sample
from iggtools.params.schemas import BLAST_M8_SCHEMA, MARKER_INFO_SCHEMA, species_profile_schema, format_data
from iggtools.params.inputs import marker_genes_hmm_cutoffs
from iggtools.params import outputs

DEFAULT_WORD_SIZE = 28
DEFAULT_ALN_COV = 0.75
DEFAULT_ALN_MAPID = 94.0
DECIMALS = ".6f"


def register_args(main_func):
    subparser = add_subcommand('midas_run_species', main_func, help='estimate species abundance profile for given sample')
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

    subparser.add_argument('--word_size',
                           dest='word_size',
                           default=DEFAULT_WORD_SIZE,
                           type=int,
                           metavar="INT",
                           help=f"Word size for BLAST search ({DEFAULT_WORD_SIZE}).  Use word sizes > 16 for greatest efficiency.")
    subparser.add_argument('--aln_mapid',
                           dest='aln_mapid',
                           type=float,
                           metavar="FLOAT",
                           help=f"Discard reads with alignment identity < ALN_MAPID.  Values between 0-100 accepted.  By default gene-specific species-level cutoffs are used, as specifeid in {marker_genes_hmm_cutoffs}")
    subparser.add_argument('--aln_cov',
                           dest='aln_cov',
                           default=DEFAULT_ALN_COV,
                           type=float,
                           metavar="FLOAT",
                           help=f"Discard reads with alignment coverage < ALN_COV ({DEFAULT_ALN_COV}).  Values between 0-1 accepted.")
    subparser.add_argument('--max_reads',
                           dest='max_reads',
                           type=int,
                           metavar="INT",
                           help=f"Number of reads to use from input file(s).  (All)")
    subparser.add_argument('--local_dbsdir',
                           dest='local_dbsdir',
                           type=str,
                           metavar="STR",
                           help=f"Provide local path of the dbs instead of sample-specific")
    return main_func


def readfq(fp):
    """ https://github.com/lh3/readfq/blob/master/readfq.py
        A generator function for parsing fasta/fastq records """
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last:
            break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield (name, ''.join(seqs), None)  # yield a fasta record
            if not last:
                break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield (name, seq, ''.join(seqs))  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield (name, seq, None)  # yield a fasta record instead
                break


def parse_reads(filename, max_reads=None):
    if not filename:
        return
    read_count_filter = None
    if max_reads != None:
        read_count_filter = f"head -n {4 * max_reads}"
    read_count = 0
    with InputStream(filename, read_count_filter) as fp:
        for name, seq, _ in readfq(fp):
            read_count += 1
            new_name = construct_queryid(name, len(seq))  # We need to encode the length in the query id to be able to recover it from hs-blastn output
            yield (new_name, seq)
        if read_count_filter:
            fp.ignore_errors()
    tsprint(f"parse_reads:: parsed {read_count} reads from {filename}")


def map_reads_hsblast(m8_file, r1, r2, word_size, markers_db, max_reads):
    assert os.path.exists(os.path.dirname(m8_file)), f"{m8_file} doesn't exit ."

    blast_command = f"hs-blastn align -word_size {word_size} -query /dev/stdin -db {markers_db} -outfmt 6 -num_threads {num_physical_cores} -evalue 1e-3"
    with OutputStream(m8_file, through=blast_command) as blast_input:
        for qid, seq in chain(parse_reads(r1, max_reads), parse_reads(r2, max_reads)):
            blast_input.write(">" + qid + "\n" + seq + "\n")


def deconstruct_queryid(rid):
    qid, qlen = rid.rsplit('_', 1)
    return qid, int(qlen)


def construct_queryid(qid, qlen):
    return f"{qid}_{qlen}"


def query_coverage(aln):
    """ Compute alignment coverage of query """
    _, qlen = deconstruct_queryid(aln['query'])
    return aln['aln'] / qlen


def find_best_hits(marker_info, m8_file, marker_cutoffs, args):
    """ Find top scoring alignment for each read """

    best_hits = {}
    i = 0
    with InputStream(m8_file) as m8_stream:
        for aln in select_from_tsv(m8_stream, schema=BLAST_M8_SCHEMA, result_structure=dict):
            i += 1
            cutoff = args.aln_mapid
            if cutoff == None:
                marker_id = marker_info[aln['target']]['marker_id'] # get gene family from marker_info
                cutoff = marker_cutoffs[marker_id]
            if aln['pid'] < cutoff: # does not meet marker cutoff
                continue
            if query_coverage(aln) < args.aln_cov: # filter local alignments
                continue
            if aln['query'] not in best_hits: # record aln
                best_hits[aln['query']] = [aln]
            elif best_hits[aln['query']][0]['score'] == aln['score']: # add aln
                best_hits[aln['query']] += [aln]
            elif best_hits[aln['query']][0]['score'] < aln['score']: # update aln
                best_hits[aln['query']] = [aln]
    tsprint(f"  total alignments: {i}")
    return list(best_hits.values())


def assign_unique(alns, species_info, marker_info):
    """ Count the number of uniquely mapped reads to each genome species """
    unique_alns = {si: [] for si in species_info}
    unique = 0
    non_unique = 0
    for aln in alns:
        if len(aln) == 1:
            unique += 1
            species_id = marker_info[aln[0]['target']]['species_id']
            unique_alns[species_id].append(aln[0])
        else:
            non_unique += 1
    tsprint(f" uniquely mapped reads: {unique}")
    tsprint(f"  ambiguously mapped reads: {non_unique}")
    return unique_alns


def assign_non_unique(alns, unique_alns, marker_info):
    """ Probabalistically assign ambiguously mapped reads """
    total_alns = unique_alns.copy()
    for aln in alns:
        if len(aln) > 1:
            species_ids = [marker_info[aln_item['target']]['species_id'] for aln_item in aln]
            counts = [len(unique_alns[sid]) for sid in species_ids]
            if sum(counts) == 0:
                species_id = random.sample(species_ids, 1)[0]
            else:
                probs = [float(count)/sum(counts) for count in counts]
                species_id = np.random.choice(species_ids, 1, p=probs)[0]
            total_alns[species_id].append(aln[species_ids.index(species_id)])
    return total_alns


def read_marker_info_repgenomes(map_file):
    with InputStream(map_file) as map_file_stream:
        return {r['gene_id']: r for r in select_from_tsv(map_file_stream, schema=MARKER_INFO_SCHEMA, result_structure=dict)}


def sum_marker_gene_lengths(marker_info):
    """ Sum up total gene length of the marker genes in each species_id's representative genome """
    total_gene_length = defaultdict(int)
    for r in marker_info.values():
        total_gene_length[r['species_id']] += int(r['gene_length'])
    return total_gene_length


def normalize_counts(species_alns, total_gene_length):
    """ Normalize counts by gene length and sum contrain """
    # norm by gene length, compute cov
    species_abundance = {}
    for species_id, alns in species_alns.items():
        # compute coverage
        if alns:
            bp = sum(aln['aln'] for aln in alns)
            cov = float(bp)/total_gene_length[species_id]
        else:
            cov = 0.0
        # TODO:  Use NamedTuple instead of dict
        species_abundance[species_id] = {'count':len(alns), 'cov':cov, 'rel_abun': 0.0}
    # compute relative abundance
    total_cov = sum(sav['cov'] for sav in species_abundance.values())
    if total_cov > 0:
        for sav in species_abundance.values():
            sav['rel_abun'] = sav['cov'] / total_cov
    tsprint(f"  total marker-gene coverage {total_cov:.3f}")
    return species_abundance


def write_abundance(species_profile_path, species_abundance):
    """ Write species results to specified output file """
    with OutputStream(species_profile_path) as outfile:
        outfile.write('\t'.join(species_profile_schema.keys()) + '\n')
        output_order = sorted(species_abundance.keys(), key=lambda sid: species_abundance[sid]['count'], reverse=True)
        for species_id in output_order:
            values = species_abundance[species_id]
            if values['count'] > 0:
                record = [species_id, values['count'], values['cov'], values['rel_abun']]
                outfile.write("\t".join(map(format_data, record)) + "\n")


def midas_run_species(args):

    try:
        sample = Sample(args.sample_name, args.midas_outdir, "species")
        sample.create_dirs(["outdir", "tempdir", "dbsdir"], args.debug)

        # Fetch db-related file either from S3 or create symlink
        dbsdir = sample.get_target_layout("dbsdir")
        if args.local_dbsdir:
            curr_dbsdir = args.local_dbsdir
            command(f"ln -s {curr_dbsdir}/* {dbsdir}")
            markers_db_files = sample.get_target_layout("marker_genes_file")
            local_toc = sample.get_target_layout("local_toc")
        else:
            markers_db_files = fetch_marker_genes(dbsdir)
            local_toc = download_reference(outputs.genomes, dbsdir)

        # Align reads to marker-genes database
        m8_file = sample.get_target_layout("species_alignments_m8")
        map_reads_hsblast(m8_file, args.r1, args.r2, args.word_size, markers_db_files[0], args.max_reads)

        with InputStream(marker_genes_hmm_cutoffs) as cutoff_params:
            marker_cutoffs = dict(select_from_tsv(cutoff_params, selected_columns={"marker_id": str, "marker_cutoff": float}))

        # Classify reads
        species_info = UHGG(local_toc).species
        marker_info = read_marker_info_repgenomes(markers_db_files[-1])
        best_hits = find_best_hits(marker_info, m8_file, marker_cutoffs, args)
        unique_alns = assign_unique(best_hits, species_info, marker_info)
        species_alns = assign_non_unique(best_hits, unique_alns, marker_info)

        # Estimate species abundance
        total_gene_length = sum_marker_gene_lengths(marker_info)
        species_abundance = normalize_counts(species_alns, total_gene_length)

        write_abundance(sample.get_target_layout("species_summary"), species_abundance)
        tsprint("Finished midas_run_species for %s" % sample.sample_name)

    except:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error. Specify --debug flag to keep.")
            sample.remove_dirs(["outdir", "tempdir", "dbsdir"])
        # TODO: find a more robust way to existing file for symlink L263
        if args.local_dbsdir:
            sample.remove_dirs(["dbsdir"])
        raise


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_run_species(args)
