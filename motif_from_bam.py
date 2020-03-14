#!/usr/bin/env python3

import pysam
import numpy

sam = pysam.AlignmentFile('sorted.allreads.bam')

motif_seqs = []
counts = numpy.zeros((4,11))



for a in sam.fetch('MT007544.1', start = 65, stop = 75):
	# ref_motif_coord = 64-a.reference_start
	# query_motif_coord = ref_motif_coord + a.query_alignment_start
	# motif_seqs.append(a.query_sequence[query_motif_coord: query_motif_coord + 11])
	
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
	if len(motif_segment) == 11:
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

