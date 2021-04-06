#!/usr/bin/env python3
import json
import os
import random
from collections import defaultdict
from itertools import chain, repeat
import numpy as np
import Bio.SeqIO

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, num_physical_cores, InputStream, OutputStream, select_from_tsv
from iggtools.models.midasdb import MIDAS_DB
from iggtools.models.sample import Sample
from iggtools.params.schemas import BLAST_M8_SCHEMA, MARKER_INFO_SCHEMA, species_profile_schema, format_data, DECIMALS6

DEFAULT_WORD_SIZE = 28
DEFAULT_ALN_COV = 0.75
DEFAULT_ALN_MAPID = 94.0


def register_args(main_func):
    subparser = add_subcommand('midas_run_species', main_func, help='Estimate species abundance profile for given sample')

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

    subparser.add_argument('--midas_db',
                           dest='midas_db',
                           type=str,
                           metavar="CHAR",
                           help=f"local MIDAS DB which mirrors the s3 IGG db")

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
                           help=f"Discard reads with alignment identity < ALN_MAPID.  Values between 0-100 accepted.  By default gene-specific species-level cutoffs are used, as specifeid in marker_genes.mapping_cutoffs.")
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
    subparser.add_argument('--num_cores',
                           dest='num_cores',
                           type=int,
                           metavar="INT",
                           default=num_physical_cores,
                           help=f"Number of physical cores to use ({num_physical_cores})")
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


def map_reads_hsblast(m8_file, r1, r2, word_size, markers_db, max_reads, num_cores):
    assert os.path.exists(os.path.dirname(m8_file)), f"{m8_file} doesn't exit."

    blast_command = f"hs-blastn align -word_size {word_size} -query /dev/stdin -db {markers_db} -outfmt 6 -num_threads {num_cores} -evalue 1e-3"
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
    tsprint(f" ambiguously mapped reads: {non_unique}")
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
    """ Only when the map file only include representative marker genes"""
    with InputStream(map_file) as map_file_stream:
        return {r['gene_id']: r for r in select_from_tsv(map_file_stream, schema=MARKER_INFO_SCHEMA, result_structure=dict)}


def read_marker_info(fasta_file, map_file):
    """ Extract gene_id_is_marker - marker_id mapping from fasta and map files """
    # Read the gene_is_marker_id from phyeco.fa file
    info = {}
    with InputStream(fasta_file) as file:
        for rec in Bio.SeqIO.parse(file, 'fasta'):
            info[rec.id] = None

    # Extract the corresponding marker_id, genome_id, and species_id from the phyeco.map file
    with InputStream(map_file) as map_file_stream:
        for r in select_from_tsv(map_file_stream, schema=MARKER_INFO_SCHEMA, result_structure=dict):
            if r['gene_id'] in info:
                info[r["gene_id"]] = r
    return info


def sum_marker_gene_lengths(marker_info):
    """ Compute the total marker gene lengths per species """
    total_gene_length = defaultdict(int)
    list_of_marker_genes = defaultdict(list)
    for r in marker_info.values():
        total_gene_length[r['species_id']] += int(r['gene_length'])
        list_of_marker_genes[r['species_id']].append(r['gene_id'])
    return total_gene_length, list_of_marker_genes


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
                #outfile.write("\t".join(map(format_data, record)) + "\n")
                outfile.write("\t".join(map(format_data, record, repeat(DECIMALS6, len(record)))) + "\n")


def midas_run_species(args):

    try:
        sample = Sample(args.sample_name, args.midas_outdir, "species")
        sample.create_dirs(["outdir", "tempdir"], args.debug)

        tsprint(f"CZ::fetch_iggdb_files::start")
        midas_db = MIDAS_DB(args.midas_db if args.midas_db else sample.get_target_layout("midas_db_dir"), 1)
        marker_db_files = midas_db.fetch_files("marker_db")
        marker_db_hmm_cutoffs = midas_db.fetch_files("marker_db_hmm_cutoffs")
        tsprint(f"CZ::fetch_iggdb_files::finish")

        with InputStream(marker_db_hmm_cutoffs) as cutoff_params:
            marker_cutoffs = dict(select_from_tsv(cutoff_params, selected_columns={"marker_id": str, "marker_cutoff": float}))

        # Align reads to marker-genes database
        tsprint(f"CZ::map_reads_hsblast::start")
        m8_file = sample.get_target_layout("species_alignments_m8")
        map_reads_hsblast(m8_file, args.r1, args.r2, args.word_size, marker_db_files["fa"], args.max_reads, args.num_cores)
        tsprint(f"CZ::map_reads_hsblast::finish")

        # Classify reads
        species_info = midas_db.uhgg.species
        #marker_info = read_marker_info_repgenomes(marker_db_files["map"])
        marker_info = read_marker_info(marker_db_files["fa"], marker_db_files["map"])
        tsprint(f"CZ::find_best_hits::start")
        best_hits = find_best_hits(marker_info, m8_file, marker_cutoffs, args)
        tsprint(f"CZ::find_best_hits::finish")

        tsprint(f"CZ::assign_unique::start")
        unique_alns = assign_unique(best_hits, species_info, marker_info)
        tsprint(f"CZ::assign_unique::finish")

        tsprint(f"CZ::assign_non_unique::start")
        species_alns = assign_non_unique(best_hits, unique_alns, marker_info)
        tsprint(f"CZ::assign_non_unique::finish")

        # Estimate species abundance
        tsprint(f"CZ::normalize_counts::start")
        total_gene_length, list_of_marker_genes = sum_marker_gene_lengths(marker_info)
        species_abundance = normalize_counts(species_alns, total_gene_length)
        tsprint(f"CZ::normalize_counts::finish")

        outfile = sample.get_target_layout("species_summary")
        write_abundance(outfile, species_abundance)

        # TODO: remove the following once finish debug
        spids = [k for k, v in species_abundance.items() if v['count'] > 0 ]
        with OutputStream(os.path.dirname(outfile) + '/marker_genes') as stream:
            for k, d, in list_of_marker_genes.items():
                if k in spids:
                    for v in d:
                        stream.write("\t".join([k, v]) + "\n")

        with OutputStream(os.path.dirname(outfile) + '/marker_length') as stream:
            for k, v, in total_gene_length.items():
                if k in spids:
                    stream.write("\t".join([k, str(v)]) + "\n")

    except Exception as error:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error. Specify --debug flag to keep.")
            sample.remove_dirs(["outdir", "tempdir"])
        raise error


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_run_species(args)
