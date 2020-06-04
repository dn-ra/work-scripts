#!/bin/bash

#conda activate ~/nanoraw_env

#inputs?
#	raw_reads_basedir $1
#	reference genome $2
#	bascalled reads $3
#	statistics-file-basename $4

#store parameters
echo ${@} > $4_params.txt

echo 'minimap'
/sw/minimap2/minimap2-2.11_x64-linux/minimap2 $2 $3 -x map-ont --secondary=no | cut -f 1 | seqtk subseq $3 - > "$4"_viralreads.fastq

echo 'rasusa'
rasusa_reads=rasusa100_"$4".fastq
rasusa --coverage 100 --genome-size 30kb --input "$4"_viralreads.fastq > $rasusa_reads

echo 'fast5 fetcher'
raw_dir=rasusa100_"$4"_fast5
mkdir $raw_dir

python3 ~/TEST_squiggle/SquiggleKit/fast5_fetcher_multi.py -q $rasusa_reads  \
-s $(dirname $3)/sequencing_summary.txt -o $raw_dir -m $1

echo 'tombo preprocess'
tombo preprocess annotate_raw_with_fastqs --overwrite --fast5-basedir $raw_dir --fastq-filenames $rasusa_reads
echo 'echo resquiggle'
tombo resquiggle --rna $raw_dir $2
cho 'tombo detect modifications'
tombo detect_modifications alternative_model --rna --fast5-basedirs $raw_dir --statistics-file-basename $4 --alternate-bases 5mC
echo 'tombo text output'
tombo text_output browser_files --browser-file-basename $4 --file-types dampened_fraction  \
--statistics-filename $4.5mC.tombo.stats --fast5-basedirs $raw_dir
