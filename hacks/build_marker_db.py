#!/usr/bin/env python3
# Find Phyeco marker genes for each genome
# Chunyu Zhao

import os
import sys
import Bio.SeqIO
import sys
import subprocess
import shutil
import gzip
import argparse
from midas import utility


def parse_fasta(p_in):
 """" Lookup of seq_id to sequence for PATRIC genes """
 seqs = {}
 infile = utility.iopen(p_in)
 for r in Bio.SeqIO.parse(infile, "fasta"):
  seqs[r.id] = str(r.seq).upper()
 infile.close()
 return seqs

def parse_hmmsearch(p_in):
 """ Parse HMMER domblout files. Return data-type formatted dictionary """
 f_in = utility.iopen(p_in)
 for line in f_in:
  if line[0] == '#': continue
  x = line.rstrip().split()
  query = x[0]
  target = x[3]
  evalue = float(x[12])
  qcov = (int(x[20]) - int(x[19]) + 1)/float(x[2])
  tcov = (int(x[16]) - int(x[15]) + 1)/float(x[5])
  yield {'query':query, 'target':target, 'evalue':evalue, 'qcov':qcov, 'tcov':tcov, 'qlen':int(x[2]), 'tlen':int(x[5])}

def find_hits(inpath, max_evalue, min_cov):
 hits = {}
 for r in parse_hmmsearch(inpath):
  if r['evalue'] > max_evalue:
   continue
  elif min(r['qcov'], r['tcov']) < min_cov:
   continue
  if r['target'] not in hits:
   hits[r['target']] = r
  elif r['evalue'] < hits[r['target']]['evalue']:
   hits[r['target']] = r
 return list(hits.values())


def main():
 usage ="""
  USAGE:
  python -genome genome -hmm hmmfile --hmm-cutoff phyeco.maping_cutoffs -mapfile mapfile -max-evalue 1e-5 -min-cov 0.00
  """

 p = argparse.ArgumentParser(
  formatter_class = argparse.RawTextHelpFormatter,
  add_help = True,
  usage = argparse.SUPPRESS,
  description = """Description:
  This script take the hmmsearch file as input, parse it based on gene-specific
  cutoffs, and build the marker genes
  """)
 ### Input
 p.add_argument('-genome', dest='genome_id', action='store', type=str, required = True,
                help='genome name.')
 p.add_argument('-genes', dest='genes_fna', action='store', type=str, required = True,
                help='ffn file for the genes annotated by Prodigal.')
 p.add_argument('-hmm', dest='hmm_file', action='store', type=str, required = True,
                help='hmm search file.')
 p.add_argument('-hmm-cutoff', dest='hmm_cutoff', action='store', type=str, default = "phyeco.mapping_cutoffs",
                help='HMM gene specific cutoffs file.')
 p.add_argument('-mapfile', dest='mapfile', action='store', type=str, default = "mapfile",
                help='genome to representative genomes mapping file.')
 p.add_argument('-max-evalue', dest='max_evalue', action='store', default=1e-5, type=float,
                help='maximum evalue for hmm search.')
 p.add_argument('-min-cov', dest='min_cov', action='store', default=0.00, type=float,
                help='minimum coverage for hmm search.')
 ### Output
 p.add_argument('-output-fasta', dest='output_fa', action='store', type=str,
                help='output fasta for the marker genes.')
 p.add_argument('-output-map', dest='output_map', action='store', type=str,
                help='output map file for the marker genes.')


 args = vars(p.parse_args())

 genome_id = args['genome_id']
 genes = args['genes_fna']
 hmm_file = args['hmm_file']
 hmm_cutoff = args['hmm_cutoff']
 mapfile = args['mapfile']
 max_evalue = args['max_evalue']
 min_cov = args['min_cov']
 output_fa = args['output_fa']
 output_map = args['output_map']


 awk_command = "awk -v pat=%s '$1 == pat {print}' %s" % (genome_id, mapfile)
 line = subprocess.getoutput(awk_command)
 species_id = line.split('\t')[1]
 is_rep = int(line.split('\t')[2])

 ffn = parse_fasta(genes)

 with open(output_fa, "w") as ofasta, open(output_map, "w") as oinfo:
  for h in find_hits(hmm_file, max_evalue, min_cov):
   gene = ffn[h['query']].upper()
   info = [species_id, genome_id, h['query'], len(gene), h['target']]
   oinfo.write('\t'.join([str(_) for _ in info]) + '\n')
   ofasta.write('>%s\n%s\n' % (h['query'], gene))



if __name__ == "__main__":
 main()
