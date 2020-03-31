#!/usr/bin/env python3

import pysam
import numpy
import argparse

#parse arguments
parser = argparse.ArgumentParser(description = 'Extract alignmened segments of reads from bam file')
required_args = parser.add_argument_group('Positional arguments')
required_args.add_argument('BAM_in', type = str, help = 'Your minibam file containing reads mapping to a certain region')


settings_args = parser.add_argument_group('Settings')
settings_args.add_argument('--jaspar', action = store_true, help = 'Do you want a fasta file of the selected region? Or a count matrix in jaspar format?')
settings_args.add_argument('--fasta', action = store_true, )
settings_args.add_argument('--coords', type = integer, nargs = 2, help = 'start and stop co-ordinates for your region of interest')
settings_args.add_argument('--chromosome', type = str, help = 'The chromosome/contig you are using as a reference', default = None)


args = parser.parse_args()

print(args)

if args.fasta == False and args.jaspar == False:
	Error('Need at least one jaspar or fasta to be selected')
	exit(1)

#-------------main------------

sam = pysam.AlignmentFile(args.BAM_in)

motif_seqs = []
motif_length = args.coords[1] - args.coords[0]
counts = numpy.zeros((4,motif_length))



for a in sam.fetch(contig = args.chromosome, start = args.coords[0], stop = args.coords[1]):
	# ref_motif_coord = 64-a.reference_start
	# query_motif_coord = ref_motif_coord + a.query_alignment_start
	# motif_seqs.append(a.query_sequence[query_motif_coord: query_motif_coord + 11])

#if args.
	pairs = a.get_aligned_pairs()
	pairs_iter = iter(pairs)
	motif_segment = ''

	for p in pairs:
		try:
			if p[1] > 75:
				break
		except TypeError:
			continue
		if p[1] >= 64 and p[1] <= 74:
            		if p[0] != None:
                		motif_segment+=a.query_sequence[p[0]]
            		else:
                		motif_segment+='-'
	if len(motif_segment) == motif_length:
	#	print(motif_segment)
		for i, sym in enumerate(motif_segment):
			if sym == 'A':
				counts[0,i] +=1
			if sym == 'C':
				counts[1,i]+=1
			if sym == 'G':
				counts[2,i]+=1
			if sym == 'T':
				counts[3,i]+=1

with open('out.jaspar', 'w') as f:
	for line in counts:
		f.write(str([int(i) for i in line])+'\n')
