#!/usr/bin/env python3

import argparse
import sys
from Bio import SeqIO
import re
import csv

parser = argparse.ArgumentParser()
parser.add_argument('--TRS_short', metavar='STRING',  type=str,help='Conserved Transcription Regulating Sequence to find. Default= ACGAAC ', default='ACGAAC', required=False)
parser.add_argument('--virus_genome', metavar='FASTA',type=str, help='Location of viral genome fasta file', required=True)
args = parser.parse_args()

genome = SeqIO.read(args.virus_genome, 'fasta')


loc = genome.seq.find(args.TRS_short, end=99)
if loc == -1:
    sys.exit('Error: Failed to find given TRS sequence in first 100bp of genome')

#location of TRS-L
start, stop = [loc+1, loc+len(args.TRS_short)] #correct for 0-based index

#exact matching for short TRS-B. write straight to file

with open('TRS_short.coords.csv', 'w') as f:
    f_write = csv.writer(f, delimiter=',', quoting=csv.QUOTE_NONE)
    f_write.writerow(['Start', 'End', 'Sequence'])
    for t in re.finditer(args.TRS_short, str(genome.seq)):
        f_write.writerow([t.span(0)[0]+1,t.span(0)[1], args.TRS_short])

