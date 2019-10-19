#!/usr/bin/env python3
# Find Phyeco marker genes for each genome
# Chunyu Zhao

import os
import sys
import Bio.SeqIO
import subprocess
import shutil
import gzip
import argparse


def print_clean_genes(gene_file, output_dir):
    gene_name = os.path.basename(gene_file)
    genome_name = gene_name.split('.fna')[0]
    output_gene = "%s/%s.genes.fna" % (output_dir, genome_name)
    output_info = "%s/%s.genes.lens" % (output_dir, genome_name)

    genes = {}
    info = {}
    for rec in Bio.SeqIO.parse(gene_file, 'fasta'):
        if str(rec.seq) == '' or str(rec.id) in ['', '|']:
            continue
        else:
            genes[rec.id] = str(rec.seq).upper()
            info[rec.id] = len(rec.seq) #gene_length

    with open(output_gene, 'w') as file:
        for geneid, geneseq in genes.items():
            file.write('>%s\n%s\n' % (geneid, geneseq))
    with open(output_info, 'w') as file:
        for geneid, genelen in info.items():
            file.write('%s\t%s\t%s\n' % (geneid, genome_name, genelen))

if __name__ == "__main__":
    p = argparse.ArgumentParser(
    formatter_class = argparse.RawTextHelpFormatter,
    add_help = True, usage = argparse.SUPPRESS,
    description = """Description:
    This script take the hmmsearch file as input, parse it based on gene-specific
    cutoffs, and build the marker genes
    """)

    p.add_argument('-gene-file', dest='gene_file', action='store', type=str, required = True,
                help='Prodigal annotated gene fna file.')
    p.add_argument('-output-dir', dest='output_dir', action='store', type=str,
                   help='output directory for single genome info.')

    args = vars(p.parse_args())
    gene_file = args['gene_file']
    output_dir = args['output_dir']
    print_clean_genes(gene_file, output_dir)
