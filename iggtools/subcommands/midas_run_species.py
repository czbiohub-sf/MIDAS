import json
import random
from collections import defaultdict
import numpy as np

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, num_physical_cores, command, InputStream, OutputStream, parse_table_as_rowdicts, parse_table, multithreading_map, download_reference, TimedSection
from iggtools import params


DEFAULT_WORD_SIZE = 28
DEFAULT_ALN_COV = 0.75


def register_args(main_func):
    subparser = add_subcommand('midas_run_species', main_func, help='compute species abundance profile for given sample(s)')
    subparser.add_argument('outdir',
                           type=str,
                           help="""Path to directory to store results.  Name should correspond to unique sample identifier.""")
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
    subparser.add_argument('--mapid',
                           dest='mapid',
                           type=float,
                           metavar="FLOAT",
                           help=f"Discard reads with alignment identity < MAPID.  Values between 0-100 accepted.  By default gene-specific species-level cutoffs are used, as specifeid in {params.inputs.marker_genes_hmm_cutoffs}")
    subparser.add_argument('--aln_cov',
                           dest='aln_cov',
                           default=DEFAULT_ALN_COV,
                           type=float,
                           metavar="FLOAT",
                           help=f"Discard reads with alignment coverage < ALN_COV ({DEFAULT_ALN_COV}).  Values between 0-1 accepted.")
    subparser.add_argument('--read_length',
                           dest='read_length',
                           type=int,
                           metavar="INT",
                           help=f"Trim reads to READ_LENGTH and discard reads with length < READ_LENGTH.  By default, reads are not trimmed or filtered")
    return main_func


def map_reads_hsblast(tempdir, r1, r2, word_size, markers_db):
    m8_file = f"{tempdir}/alignments.m8"
    cat_cmd = f"cat {r1}"
    if r2:
        cat_cmd += f" {r2}"
    # TODO: Strip irrelvant parts of read ID, see stream_seq in MIDAS
    command(f"{cat_cmd} | hs-blastn align -word_size {word_size} -query /dev/stdin -db {markers_db} -outfmt 6 -num_threads {num_physical_cores} -out {m8_file} -evalue 1e-3")
    return m8_file


def parse_blast(inpath):
    """ Yield formatted record from BLAST m8 file """
    # TODO: Bring schema parser from idseq-dag.
    formats = [str, str, float, int, float, float, float, float, float, float, float, float]
    fields = ['query', 'target', 'pid', 'aln', 'mis', 'gaps', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'score']
    for line in open(inpath):
        values = line.rstrip().split()
        yield dict((field, format(value)) for field, format, value in zip(fields, formats, values))


def query_coverage(aln):
    """ Compute alignment coverage of query """
    qlen = aln['query'].split('_')[-1] # get qlen from sequence header
    return float(aln['aln']) / int(qlen)


def find_best_hits(args, marker_info, m8_file, marker_cutoffs):
    """ Find top scoring alignment for each read """
    best_hits = {}
    i = 0
    for aln in parse_blast(m8_file):
        i += 1
        cutoff = args.mapid
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
            #species_id = aln[0]['target'].split('_')[0]
            species_id = marker_info[aln[0]['target']]['species_id']
            unique_alns[species_id].append(aln[0])
        else:
            non_unique += 1
    tsprint(f"  uniquely mapped reads: {unique}")
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
    columns = ["species_id", "genome_id", "gene_id", "gene_length", "marker_id"]
    with InputStream(map_file) as impf:
        return {r['gene_id']: r for r in parse_table_as_rowdicts(impf, columns, columns)}


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


def write_abundance(outdir, species_abundance):
    """ Write species results to specified output file """
    outpath = f"{outdir}/species/species_profile.txt"
    with OutputStream(outpath) as outfile:
        fields = ['species_id', 'count_reads', 'coverage', 'relative_abundance']
        outfile.write('\t'.join(fields) + '\n')
        output_order = sorted(species_abundance.keys(), key=lambda sid: species_abundance[sid]['count'], reverse=True)
        for species_id in output_order:
            values = species_abundance[species_id]
            record = [species_id, values['count'], values['cov'], values['rel_abun']]
            outfile.write('\t'.join(str(x) for x in record) + '\n')


def midas_run_species(args):

    tempdir = f"{args.outdir}/species/temp/"

    command(f"rm -rf {tempdir}")
    command(f"mkdir -p {tempdir}")

    markers_db_files = multithreading_map(download_reference, ("s3://microbiome-igg/2.0/marker_genes/phyeco/phyeco.fa{ext}.lz4" for ext in ["", ".bwt", ".header", ".sa", ".sequence"]))

    local_toc = download_reference(params.outputs.genomes)
    species_info, _, _ = read_toc(local_toc)

    marker_info = read_marker_info_repgenomes(...)

    with TimedSection("aligning reads to marker-genes database"):
        m8_file = map_reads_hsblast(tempdir, args.r1, args.r2, args.word_size, markers_db_files[0])

    with InputStream(params.inputs.marker_genes_hmm_cutoffs) as cutoff_params:
        marker_cutoffs = dict(parse_table(cutoff_params, ["marker_id", "marker_cutoff"]))

    with TimedSection("classifying reads"):
        best_hits = find_best_hits(args, marker_info, m8_file, marker_cutoffs)
        unique_alns = assign_unique(best_hits, species_info, marker_info)
        species_alns = assign_non_unique(best_hits, unique_alns, marker_info)

    with TimedSection("estimating species abundance"):
        total_gene_length = sum_marker_gene_lengths(marker_info)
        species_abundance = normalize_counts(species_alns, total_gene_length)

    write_abundance(args.outdir, species_abundance)


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args\n{json.dumps(vars(args), indent=4)}")
    midas_run_species(args)
