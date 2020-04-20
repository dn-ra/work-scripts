#!/usr/bin/env python3

import argparse
import sys
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--TRS_seq', metavar='STRING',  type=str,help='Conserved Transcription Regulating Sequence to find. Default= ACGAAC ', default='ACGAAC', required=False)
parser.add_argument('--virus_genome', metavar='FASTA',type=str, help='Location of viral genome fasta file', required=True)
args = parser.parse_args()

genome = SeqIO.read(args.virus_genome, 'fasta')[0:99]

loc = genome.seq.find(args.TRS_seq, end=99)
if loc == -1:
    sys.exit('Error: Failed to find given TRS sequence in first 100bp of genome')
else:
    start, stop = [loc+1, loc+len(args.TRS_seq)] #correct for 0-based index
  

        
