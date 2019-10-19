#!/usr/bin/env python3
# Chunyu Zhao

import os
import sys
import Bio.SeqIO
import argparse
from midas import utility
from collections import defaultdict

def parse_uclust(inpath):
    """ Yield formated records from UCLUST outputfile """
    fields = ['type', 'cluster_id', 'size', 'pid', 'strand', 'skip1', 'skip2', 'skip3', 'gene_id', 'centroid_id']
    with utility.iopen(inpath) as infile:
        for index, line in enumerate(infile):
            values = line.rstrip('\n').split('\t')
            record = dict([(f, v) for f, v in zip(fields, values)])
            yield record

def store_gene_info(uclust_dir, output_file):

    genes = defaultdict(lambda: defaultdict(int))
    stats = {}

    """ store_gene_info """
    for pid in [99, 95, 90, 85, 80, 75]:
        stats['centroids_%s' % pid] = 0
        uclust_file = '%s/uclust.%s.txt' % (uclust_dir, pid)
        for r in parse_uclust(uclust_file):
            if r['type'] == 'H': ## non-centroids
                genes[r['gene_id']]["cluster_id_%s" % str(pid)] = r['cluster_id']
                genes[r['gene_id']]["centroid_id_%s" % str(pid)] = r['centroid_id']
            elif r['type'] == 'S': ## centroids
                genes[r['gene_id']]["cluster_id_%s" % str(pid)] = r['cluster_id']
                genes[r['gene_id']]["centroid_id_%s" % str(pid)] = r['gene_id']
                stats['centroids_%s' % pid] += 1
            else:
                continue

    """ store_cluster_membership: Map gene to 99% ID centroids at each clustering %ID cutoff """
    for gene in genes.values():
        gene["centroid_99"] = gene["centroid_id_99"]
        gene["centroid_95"] = genes[ gene["centroid_99"] ]["centroid_id_95"]
        gene["centroid_90"] = genes[ gene["centroid_99"] ]["centroid_id_90"]
        gene["centroid_85"] = genes[ gene["centroid_99"] ]["centroid_id_85"]
        gene["centroid_80"] = genes[ gene["centroid_99"] ]["centroid_id_80"]
        gene["centroid_75"] = genes[ gene["centroid_99"] ]["centroid_id_75"]

    """ write_gene_info """
    file = utility.iopen(output_file, 'w')
    header = ['gene_id', 'centroid_99', 'centroid_95', 'centroid_90', 'centroid_85', 'centroid_80', 'centroid_75']
    file.write('\t'.join(header)+'\n')
    for gene_id in sorted(genes.keys()):
        g = genes[gene_id]
        values = [gene_id, g["centroid_99"], g["centroid_95"], g["centroid_90"], g["centroid_85"], g["centroid_80"], g["centroid_75"]]
        file.write('\t'.join(str(_) for _ in values) + '\n')
    file.close()


if __name__ == "__main__":
    p = argparse.ArgumentParser(
    formatter_class = argparse.RawTextHelpFormatter,
    add_help = True, usage = argparse.SUPPRESS,
    description = """Description:
    This script take the hmmsearch file as input, parse it based on gene-specific
    cutoffs, and build the marker genes
    """)

    p.add_argument('-uclust-dir', dest='uclust_dir', action='store', type=str, required = True,
                help='genome name.')
    p.add_argument('-output-file', dest='output_file', action='store', type=str,
                   help='output fasta for the marker genes.')

    args = vars(p.parse_args())

    uclust_dir = args['uclust_dir']
    output_file = args['output_file']
    store_gene_info(uclust_dir, output_file)
