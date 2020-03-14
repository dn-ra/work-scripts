#!/bin/env python3

import pysam

sam = pysam.AlignmentFile('sorted.allreads.bam')

motif_seqs = []


for a in sam.fetch('MT007544.1', start = 65, stop = 75):
	# ref_motif_coord = 64-a.reference_start
	# query_motif_coord = ref_motif_coord + a.query_alignment_start
	# motif_seqs.append(a.query_sequence[query_motif_coord: query_motif_coord + 11])
	
	pairs = a.get_aligned_pairs()
	pairs_iter = iter(pairs)
	motif_segment = []
	
	for p in pairs:
		motif_segment = ''
		if p[1] > 75:
			break
		elif p[1] >= 64 && p[1] <= 74:
		if p[0] != None:
			motif_segment+=p[0]
		else:
			motif_segment+='-'
	motif_seqs.append(motif_segment)
		
	
	