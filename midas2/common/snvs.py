#!/usr/bin/env python3
from operator import itemgetter
import numpy as np # pylint: disable=no-name-in-module
from midas2.params.schemas import format_data


def query_overlap_qualities(f, r):
    # "The higher quality base is used and the lower-quality base is set to BQ=0."
    for i, _ in enumerate(zip(f, r)):
        (x, y) = _
        if x >= y:
            r[i] = 0
        else:
            f[i] = 0
    return (f, r)


def is_allele(readcount, site_depth, snp_maf_cutoff, allele_depth_cutoff=2):
    # keep alleles passing (1) min allele read count (2) the min allele frequency
    return readcount >= allele_depth_cutoff and readcount / site_depth >= snp_maf_cutoff


def call_alleles(tuple_of_alleles, site_depth, snp_maf_cutoff):
    """ Call SNPs and compute allele frequencies """

    # Only when you have seen all the revelant samples, you can call SNPs
    alleles_above_cutoff = tuple(al for al in tuple_of_alleles if is_allele(al[1], site_depth, snp_maf_cutoff))

    # classify SNPs type
    number_alleles = len(alleles_above_cutoff)
    if number_alleles == 0:
        return (None, None, None, 0)

    snp_type = ["mono", "bi", "tri", "quad"][number_alleles - 1]

    # In the event of a tie -- biallelic site with 50/50 freq split -- the allele declared major is
    # the one that comes later in the "ACGT" lexicographic order.
    alleles_above_cutoff = sorted(alleles_above_cutoff, key=itemgetter(1), reverse=True)[:2]
    major_allele = alleles_above_cutoff[0][0]
    minor_allele = alleles_above_cutoff[-1][0] # for fixed sites, same as major allele

    return (major_allele, minor_allele, snp_type, number_alleles)


def reference_overlap(p, q):
    return max(0.0, min(p[1], q[1]) - max(p[0], q[0]) + 1)


def hamming_distance(str1, str2):
    assert len(str1) == len(str2), f"Two input strings for hamming_distance are different length."
    hd = 0
    for i in range(len(str1)):
        if str1[i] != str2[i] or str1[i] not in 'ACGT':
            hd += 1
    return hd


def position_within_overlap(pos, strand, boundary):
    """ Check if given position is within the overlap range """
    if strand == "fwd" and pos is not None and pos >= boundary:
        return True
    if strand == "rev" and pos is not None and pos <= boundary:
        return True
    return False


def update_overlap(reads_overlap, aln):
    """ The actual overlap should substract the number of gaps in the forward read """
    aligned_pos = aln.get_aligned_pairs()
    ngaps = 0
    for i in range(0, len(aligned_pos)):
        if aligned_pos[i][1] is not None and aligned_pos[i][1] >= aln.reference_end - reads_overlap and aligned_pos[i][0] is None:
            ngaps += 1
    return reads_overlap - ngaps


def mismatches_within_overlaps(aln, reads_overlap, strand):
    """ For given alignment, compute NM within and outside overlap with paired read """

    # reference sequence that is covered by reads alignment
    ref_seq = aln.get_reference_sequence()
    # aligned portion of the read
    qry_seq = aln.query_alignment_sequence
    # a list of aligned read (query) and reference positions
    aligned_pos = aln.get_aligned_pairs()

    ngaps_qi = 0
    ngaps_ri = 0
    ngaps_qo = 0
    ngaps_ro = 0
    ro = []
    qo = []
    ri = []
    qi = []

    for i in range(0, len(aligned_pos)):

        boundary = aln.query_alignment_end - reads_overlap if strand == "fwd" else aln.query_alignment_start + reads_overlap - 1

        if position_within_overlap(aligned_pos[i][0], strand, boundary):
            # The aligned position is witin the overlap
            if aligned_pos[i][0] is None:
                qi.append("-")
                ngaps_qi += 1
            else:
                qi.append(qry_seq[aligned_pos[i][0] - aln.query_alignment_start])

            if aligned_pos[i][1] is None:
                ri.append("-")
                ngaps_ri += 1
            else:
                ri.append(ref_seq[aligned_pos[i][1] - aln.reference_start])
        else:
            # The aligned position is outside the overlap
            if aligned_pos[i][0] is None:
                qo.append("-")
                ngaps_qo += 1
            else:
                qo.append(qry_seq[aligned_pos[i][0] - aln.query_alignment_start])

            if aligned_pos[i][1] is None:
                ro.append("-")
                ngaps_ro += 1
            else:
                ro.append(ref_seq[aligned_pos[i][1] - aln.reference_start])

    ro = "".join(ro)
    qo = "".join(qo)
    nm_out = hamming_distance(ro.upper(), qo.upper())

    ri = "".join(ri)
    qi = "".join(qi)
    nm_in = hamming_distance(ri.upper(), qi.upper())

    alned_no_gaps = len((ri + ro).replace("-", ""))

    ngaps_q = ngaps_qi + ngaps_qo
    ngaps_r = ngaps_ri + ngaps_ro


    row = ["func::mismatches_within_overlaps", aln.query_alignment_length, alned_no_gaps, ngaps_r, ngaps_q, aln.query_name,
           aln.reference_name, aln.reference_start, aln.reference_end, aln.query_length, aln.get_aligned_pairs()]
    if ngaps_qi > 0:
        print("\t".join(map(str, row)))

    #assert ngaps_qi == 0
    #assert aln.query_alignment_length == len((qi + qo).replace("-", ""))
    #assert aln.query_alignment_length + ngaps_q == alned_no_gaps + ngaps_r, "\n".join(map(str, row)) + "\n"

    return (nm_out, nm_in, ngaps_ri, ngaps_ro)


def debug_overlap(alns):
    aln = alns["fwd"]
    row = [aln.reference_name, aln.reference_start, aln.reference_end,
           aln.query_name, aln.query_alignment_start, aln.query_alignment_end,
           "R1" if aln.is_read1 else "R2", "rev" if aln.is_reverse else "fwd",
           dict(aln.tags)['NM'], aln.query_alignment_length, aln.query_length]
           #reads_overlap, fragment_length, mismatches, ngaps_ri_rev, ngaps_ro_rev, ngaps_ri_fwd, ngaps_ro_fwd]
    print("+++++++++++++++++++++++++++++++++++++++++++")
    print("\t".join(map(format_data, row)))
    print(aln.get_aligned_pairs())
    aln = alns["rev"]
    row = [aln.reference_name, aln.reference_start, aln.reference_end,
           aln.query_name, aln.query_alignment_start, aln.query_alignment_end,
           "R1" if aln.is_read1 else "R2", "rev" if aln.is_reverse else "fwd",
           dict(aln.tags)['NM'], aln.query_alignment_length, aln.query_length]
           #reads_overlap, fragment_length, mismatches, ngaps_ri_rev, ngaps_ro_rev, ngaps_ri_fwd, ngaps_ro_fwd]
    print("\t".join(map(format_data, row)))
    print(aln.get_aligned_pairs())
    assert False, aln.reference_name


def _print_aln(aln):
    """ For debugging purpose, print out filter related information """
    rname = aln.reference_name
    qname = aln.query_name
    align_len = len(aln.query_alignment_sequence)
    query_len = aln.query_length
    pid = 100 * (align_len - dict(aln.tags)['NM']) / float(align_len)
    readq = np.mean(aln.query_qualities)
    mapq = aln.mapping_quality
    alncov = align_len / float(query_len)
    return [rname, qname, aln.qstart, aln.qend, pid, readq, mapq, alncov, aln.is_secondary]
