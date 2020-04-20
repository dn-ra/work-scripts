#!/usr/bin/env python3

import argparse
import sys
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--TRS_seq', type=str,help='Conserved Transcription Regulating Sequence to find', default='ACGAAC', required=False)
parser.add_argument('--virus_genome', type=str, help='Location of viral genome fasta file', required=True)
args = parser.parse_args()

lead100 = SeqIO.read(args.virus_genome, 'fasta')[0:99]

loc = lead100.seq.find(args.TRS_seq)
if loc == -1:
    sys.exit('Error: Failed to find given TRS sequence in first 100bp of genome')
else:
    start, stop = [loc+1, loc+len(args.TRS_seq)+1] #correct for 0-based index


        
