#!/usr/bin/env python3
import json
import os
import random
from collections import defaultdict
from itertools import chain, repeat
import numpy as np
import Bio.SeqIO

from midas2.common.argparser import add_subcommand
from midas2.common.utils import tsprint, num_physical_cores, InputStream, OutputStream, select_from_tsv
from midas2.models.midasdb import MIDAS_DB
from midas2.models.sample import Sample
from midas2.params.schemas import BLAST_M8_SCHEMA, MARKER_INFO_SCHEMA, species_profile_schema, species_marker_profile_schema, format_data, DECIMALS6
from midas2.params.inputs import MIDASDB_NAMES


DEFAULT_WORD_SIZE = 28
DEFAULT_ALN_COV = 0.75
DEFAULT_ALN_MAPID = 94.0
DEFAULT_MARKER_READS = 2
DEFAULT_MARKER_COVERED = 2


def register_args(main_func):
    subparser = add_subcommand('run_species', main_func, help='Estimate species abundance profile for given sample')

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

    subparser.add_argument('--midasdb_name',
                           dest='midasdb_name',
                           type=str,
                           default="uhgg",
                           choices=MIDASDB_NAMES,
                           help=f"MIDAS Database name.")
    subparser.add_argument('--midasdb_dir',
                           dest='midasdb_dir',
                           type=str,
                           default="midasdb",
                           help=f"Local MIDAS Database path mirroing S3.")

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

    subparser.add_argument('--marker_reads',
                           dest='marker_reads',
                           default=DEFAULT_MARKER_READS,
                           type=float,
                           metavar="FLOAT",
                           help=f"Only report marker covered by at least {DEFAULT_MARKER_READS} reads.")
    subparser.add_argument('--marker_covered',
                           dest='marker_covered',
                           default=DEFAULT_MARKER_COVERED,
                           type=float,
                           metavar="FLOAT",
                           help=f"Only species with more than {DEFAULT_MARKER_COVERED} markeres by at least {DEFAULT_MARKER_READS} reads.")

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


def map_reads_hsblastn(m8_file, r1, r2, word_size, markers_db, max_reads, num_cores):
    assert os.path.exists(os.path.dirname(m8_file)), f"{m8_file} doesn't exit."
    blast_command = f"hs-blastn align -outfmt 6 -num_threads {num_cores} -evalue 1e-3 -word_size {word_size} -query /dev/stdin -db {markers_db}"
    #blast_command = f"hs-blastn -outfmt 6 -num_threads {num_cores} -evalue 1e-3 /dev/stdin {markers_db}"
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


def read_markers_info(fasta_file, map_file, genes_that_are_marker_fp):
    """ Extract gene_id_is_marker - marker_id mapping from fasta and map files """
    # Read the gene_is_marker_id from phyeco.fa file
    genes_that_are_marker = []
    with InputStream(fasta_file) as file:
        for rec in Bio.SeqIO.parse(file, 'fasta'):
            genes_that_are_marker.append(rec.id)

    with OutputStream(genes_that_are_marker_fp) as stream:
        stream.write("\n".join(genes_that_are_marker) + "\n")

    # Extract the mapping of <marker_id, genome_id, species_id> from the phyeco.map file
    markers_info = defaultdict(dict)
    markers_length = defaultdict(lambda: defaultdict(int))
    markers_gene_list = defaultdict(lambda: defaultdict(list))
    filter_cmd = f"grep -Fwf {genes_that_are_marker_fp}"

    with InputStream(map_file, filter_cmd) as stream:
        for r in select_from_tsv(stream, schema=MARKER_INFO_SCHEMA, result_structure=dict):
            markers_info[r["gene_id"]] = r
            markers_length[r['species_id']][r['marker_id']] += int(r["gene_length"])
            markers_gene_list[r['species_id']][r['marker_id']].append(r["gene_id"])
        stream.ignore_errors()
    return markers_info, markers_length, markers_gene_list


def find_best_hits(m8_file, markers_info, marker_cutoffs, args):
    """ Find top scoring alignment for each read """
    best_hits = {}
    i = 0
    with InputStream(m8_file) as m8_stream:
        for aln in select_from_tsv(m8_stream, schema=BLAST_M8_SCHEMA, result_structure=dict):
            i += 1

            # Default specific marker genes sequence identity cutoff
            cutoff = marker_cutoffs[markers_info[aln['target']]['marker_id']]
            cutoff = args.aln_mapid if args.aln_mapid is not None else cutoff
            if aln['pid'] < cutoff:
                continue

            if query_coverage(aln) < args.aln_cov: # filter local alignments
                continue

            # For each read (query), find the best hits based on the reported score
            if aln['query'] not in best_hits: # record aln
                best_hits[aln['query']] = [aln]
            elif best_hits[aln['query']][0]['score'] == aln['score']: # add aln
                best_hits[aln['query']] += [aln]
            elif best_hits[aln['query']][0]['score'] < aln['score']: # update aln
                best_hits[aln['query']] = [aln]

    tsprint(f"  total alignments: {i}")
    return list(best_hits.values())


def filter_species_by_alns(alns, min_mreads=2, min_mcounts=2):
    """ Only species with more then 2 markers covered by 2 reads are reported """
    # Each marker needs to be covered by more than one reads
    filtered_alns = defaultdict(lambda: defaultdict(dict))
    filtered_covered_markers = defaultdict(list)

    for spid, aln_dict in alns.items():
        for mkid, mkdict in aln_dict.items():
            if mkdict['readcounts'] >= min_mreads:
                filtered_alns[spid][mkid] = mkdict
                filtered_covered_markers[spid].append(mkid)

    # At least two markers are covered with at least two reads
    final_covered_markers = {spid:lom for spid, lom in filtered_covered_markers.items() if len(lom) >= min_mcounts}

    final_alns = {spid:_ for spid, _ in filtered_alns.items() if spid in final_covered_markers.keys()}

    return final_alns, final_covered_markers


def assign_unique(alns, markers_info, args):
    """
    Assign uniquely mapped read to each marker gene
    final_unique_alns are indexed by <species_id, marker_id>:
        - alns: total aligned bps
        - readcounts: totall mapped reads
    """

    unique_alns = defaultdict(lambda: defaultdict(dict))
    unique_counts = defaultdict(int)
    non_unique_counts = 0

    unique_reads = defaultdict(lambda: defaultdict(list))
    for aln in alns:
        if len(aln) == 1:
            gid = aln[0]['target']
            spid = markers_info[gid]['species_id']
            mkid = markers_info[gid]['marker_id']

            unique_reads[spid][mkid].append("@"+aln[0]['query'])

            if mkid in unique_alns[spid]:
                unique_alns[spid][mkid]["alns"] += aln[0]['aln']
                unique_alns[spid][mkid]["readcounts"] += 1
            else:
                unique_alns[spid][mkid] = {"alns": aln[0]['aln'], "readcounts": 1}

            unique_counts[spid] += 1
        else:
            non_unique_counts += 1

    # At least two markers covered by at least 2 reads each
    final_unique_alns, final_covered_markers = filter_species_by_alns(unique_alns, args.marker_reads, args.marker_covered)

    tsprint(f" uniquely mapped reads: {sum(unique_counts.values())}")
    tsprint(f" ambiguously mapped reads: {non_unique_counts}")

    return final_unique_alns, final_covered_markers


def assign_non_unique(alns, unique_alns, markers_info, args):
    """ Probabilistically assign ambiguously mapped reads to markers """

    ambiguous_alns = defaultdict(lambda: defaultdict(dict))

    for loaln in alns:
        if len(loaln) > 1:
            # Special case: when the same gene was mapped twice to different regions, we only record once.
            target_dict = defaultdict(dict) # indexed by query gene
            for aln_item in loaln:
                gid = aln_item["target"]
                spid = markers_info[gid]["species_id"]
                mkid = markers_info[gid]["marker_id"]

                uniq_count = unique_alns[spid][mkid]["readcounts"] if spid in unique_alns and mkid in unique_alns[spid] else 0
                target_dict[gid] = {"species_id": spid, "marker_id": mkid, "alns": aln_item['aln'], "uniq_count": uniq_count}

            lo_uniq_counts = [_['uniq_count'] for _ in target_dict.values()]
            if sum(lo_uniq_counts) == 0:
                gene_id = random.sample(list(target_dict.keys()), 1)[0]
            else:
                probs = [float(count)/sum(lo_uniq_counts) for count in lo_uniq_counts]
                gene_id = np.random.choice(list(target_dict.keys()), 1, p=probs)[0]

            # Probabilistically assigned: gene_id
            mkid = target_dict[gene_id]['marker_id']
            spid = target_dict[gene_id]['species_id']
            aln_bps = target_dict[gene_id]['alns']
            if mkid in ambiguous_alns[spid]:
                ambiguous_alns[spid][mkid]["alns"] += aln_bps
                ambiguous_alns[spid][mkid]["readcounts"] += 1
            else:
                ambiguous_alns[spid][mkid] = {"alns": aln_bps, "readcounts": 1}

    final_ambiguous_alns, final_covered_markers = filter_species_by_alns(ambiguous_alns, args.marker_reads, args.marker_covered)

    return final_ambiguous_alns, final_covered_markers


def merge_counts(unique_alns, ambiguous_alns, unique_covered_markers, ambiguous_covered_markers, markers_length):
    """ Merge unique and ambiguous alns into full species by marker matrix """

    list_of_all_species = set(list(unique_alns.keys()) + list(ambiguous_alns.keys()))

    species_alns = defaultdict(lambda: defaultdict(dict))
    for spid in list_of_all_species:
        list_of_marker_ids = markers_length[spid].keys()
        for mkid in list_of_marker_ids:
            uniq_count = unique_alns[spid][mkid]['readcounts'] if spid in unique_alns and mkid in unique_alns[spid] else 0
            amb_count = ambiguous_alns[spid][mkid]['readcounts'] if spid in ambiguous_alns and mkid in ambiguous_alns[spid] else 0

            uniq_bps = unique_alns[spid][mkid]['alns'] if spid in unique_alns and mkid in unique_alns[spid] else 0
            amb_bps = ambiguous_alns[spid][mkid]['alns'] if spid in ambiguous_alns and mkid in ambiguous_alns[spid] else 0

            species_alns[spid][mkid] = {"unique": uniq_count, "ambiguous": amb_count, "unique_bps": uniq_bps, "ambiguous_bps": amb_bps}

    species_covered_markers = defaultdict(dict)
    for spid, marker_dict in species_alns.items():
        lom_uniq = unique_covered_markers[spid] if spid in unique_alns else []
        lom_amb = ambiguous_covered_markers[spid] if spid in ambiguous_alns else []

        uniq_covered = len(lom_uniq)
        amb_covered = len(lom_amb)
        total_covered = len(set(lom_uniq + lom_amb))

        species_covered_markers[spid] = {"total_covered":total_covered, "unique_covered":uniq_covered, "ambiguous_covered":amb_covered, "total_marker_counts": len(marker_dict.keys())}

    return species_alns, species_covered_markers


def normalize_counts(species_alns, species_covered_markers, markers_length, markers_gene_list):
    """ Normalize counts by gene length and sum contrain """

    species_abundance = defaultdict(lambda: defaultdict(int)) # indexed by <species_id>
    markers_abundance = defaultdict(lambda: defaultdict(dict)) # indexed by <species_id, marker_id>
    for spid, sp_mkdict in species_alns.items():
        lomc = []
        for mkid, mkdict in sp_mkdict.items():
            # For each marker gene, compute coverage
            mklength = markers_length[spid][mkid]
            aln_bps = mkdict["unique_bps"] + mkdict["ambiguous_bps"]
            readcounts = mkdict["unique"] + mkdict["ambiguous"]

            mkcov = aln_bps / mklength
            lomc.append(mkcov)

            # Accumulate information across markers per species
            species_abundance[spid]["read_counts"] += readcounts
            species_abundance[spid]["total_bps"] += aln_bps
            species_abundance[spid]["total_marker_length"] += mklength

            # Record information for <species_id, marker_id>
            markers_abundance[spid][mkid]["gene_id"] = ",".join(markers_gene_list[spid][mkid])
            markers_abundance[spid][mkid]["length"] = mklength
            markers_abundance[spid][mkid]["coverage"] = mkcov
            markers_abundance[spid][mkid]["total_reads"] = mkdict["unique"] + mkdict["ambiguous"]
            markers_abundance[spid][mkid]["uniq_reads"] = mkdict["unique"]
            markers_abundance[spid][mkid]["ambi_reads"] = mkdict["ambiguous"]
            markers_abundance[spid][mkid]["total_alnbps"] = mkdict["unique_bps"] + mkdict["ambiguous_bps"]
            markers_abundance[spid][mkid]["uniq_alnbps"] = mkdict["unique_bps"]
            markers_abundance[spid][mkid]["ambi_alnbps"] = mkdict["ambiguous_bps"]

        # Second round of collecting information per species
        species_abundance[spid]["coverage"] = species_abundance[spid]["total_bps"] / species_abundance[spid]["total_marker_length"]
        species_abundance[spid]["median_coverage"] = np.median(lomc)

        species_abundance[spid]["total_covered"] = species_covered_markers[spid]["total_covered"]
        species_abundance[spid]["unique_covered"] = species_covered_markers[spid]["unique_covered"]
        species_abundance[spid]["ambiguous_covered"] = species_covered_markers[spid]["ambiguous_covered"]
        species_abundance[spid]["total_marker_counts"] = species_covered_markers[spid]["total_marker_counts"]
        species_abundance[spid]["unique_fraction_covered"] = species_abundance[spid]["unique_covered"] / species_abundance[spid]["total_marker_counts"]

    # Third round: compute relative abundance based on vertical coverage: only among detected species.
    sample_coverage = sum(_['coverage'] for _ in species_abundance.values())
    for r in species_abundance.values():
        r["relative_abundance"] = r["coverage"] / sample_coverage if sample_coverage > 0 else 0.0

    tsprint(f"  total marker-gene coverage {sample_coverage:.3f}")
    return species_abundance, markers_abundance


def write_abundance(species_path, markers_path, species_abundance, markers_abundance):
    """ Write species results to specified output file """
    # Sort the species by median_coverage
    output_order = sorted(species_abundance.keys(), key=lambda sid: species_abundance[sid]['median_coverage'], reverse=True)

    with OutputStream(species_path) as outfile:
        outfile.write('\t'.join(species_profile_schema.keys()) + '\n')
        for species_id in output_order:
            r = species_abundance[species_id]
            record = [species_id, r['read_counts'], r['median_coverage'], r['coverage'], r['relative_abundance'],
                      r['total_covered'], r['unique_covered'], r['ambiguous_covered'],
                      r['total_marker_counts'], r['unique_fraction_covered'], r['total_marker_length']]
            outfile.write("\t".join(map(format_data, record, repeat(DECIMALS6, len(record)))) + "\n")

    with OutputStream(markers_path) as outfile:
        outfile.write('\t'.join(species_marker_profile_schema.keys()) + '\n')
        for species_id in output_order:
            for mid, md in markers_abundance[species_id].items():
                record = [species_id, mid, md['length'], md["gene_id"], md['total_reads'], md['total_alnbps'],
                          md['coverage'], md['uniq_reads'], md['ambi_reads'], md['uniq_alnbps'], md['ambi_alnbps']]
                outfile.write("\t".join(map(format_data, record, repeat(DECIMALS6, len(record)))) + "\n")


def run_species(args):

    try:
        sample = Sample(args.sample_name, args.midas_outdir, "species")
        sample.create_dirs(["outdir", "tempdir"], args.debug)

        tsprint(f"MIDAS2::fetch_midasdb_files::start")
        midas_db = MIDAS_DB(os.path.abspath(args.midasdb_dir), args.midasdb_name)

        marker_db_files = midas_db.fetch_files("marker_db")
        marker_db_hmm_cutoffs = midas_db.fetch_files("marker_db_hmm_cutoffs")

        tsprint(f"MIDAS2::fetch_midasdb_files::finish")

        with InputStream(marker_db_hmm_cutoffs) as cutoff_params:
            marker_cutoffs = dict(select_from_tsv(cutoff_params, selected_columns={"marker_id": str, "marker_cutoff": float}))

        # Align reads to marker-genes database
        tsprint(f"MIDAS2::map_reads_hsblastn::start")
        m8_file = sample.get_target_layout("species_alignments_m8")
        if args.debug and os.path.exists(m8_file):
            tsprint(f"Use existing {m8_file} according to --debug flag.")
        else:
            map_reads_hsblastn(m8_file, args.r1, args.r2, args.word_size, marker_db_files["fa"], args.max_reads, args.num_cores)
        tsprint(f"MIDAS2::map_reads_hsblastn::finish")

        tsprint(f"MIDAS2::read in marker information::start")
        genes_that_are_marker_fp = sample.get_target_layout("species_marker_genes")
        markers_info, markers_length, markers_gene_list = read_markers_info(marker_db_files["fa"], marker_db_files["map"], genes_that_are_marker_fp)
        tsprint(f"MIDAS2::read in marker information::finish")

        # Classify reads
        tsprint(f"MIDAS2::find_best_hits::start")
        best_hits = find_best_hits(m8_file, markers_info, marker_cutoffs, args)
        tsprint(f"MIDAS2::find_best_hits::finish")

        tsprint(f"MIDAS2::assign_unique::start")
        unique_alns, unique_covered_markers = assign_unique(best_hits, markers_info, args)
        tsprint(f"MIDAS2::assign_unique::finish")

        tsprint(f"MIDAS2::assign_non_unique::start")
        ambiguous_alns, ambiguous_covered_markers = assign_non_unique(best_hits, unique_alns, markers_info, args)
        tsprint(f"MIDAS2::assign_non_unique::finish")

        # Estimate species abundance
        tsprint(f"MIDAS2::normalize_counts::start")
        species_alns, species_covered_markers = merge_counts(unique_alns, ambiguous_alns, unique_covered_markers, ambiguous_covered_markers, markers_length)
        species_abundance, markers_abundance = normalize_counts(species_alns, species_covered_markers, markers_length, markers_gene_list)
        tsprint(f"MIDAS2::normalize_counts::finish")

        write_abundance(sample.get_target_layout("species_summary"), sample.get_target_layout("markers_summary"), species_abundance, markers_abundance)

    except Exception as error:
        if not args.debug:
            tsprint("Deleting untrustworthy outputs due to error. Specify --debug flag to keep.")
            sample.remove_dirs(["outdir", "tempdir"])
        raise error


@register_args
def main(args):
    tsprint(f"Species abundance estimation in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    run_species(args)
