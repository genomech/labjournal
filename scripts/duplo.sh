#!/bin/bash

INPUT_FOLDER="/dev/datasets/FairWind/_results/60m/PRIMARY_ANALYSIS_13D/bam_sorted"

for name in 'sample-1-1'
do

mkdir $INPUT_FOLDER/"$name"_temp;

PicardCommandLine MarkDuplicates REMOVE_DUPLICATES=true DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES \
	M=$INPUT_FOLDER/"$name"_temp/1.txt \
	I=$INPUT_FOLDER/"$name"_60M_sorted.bam \
	O=$INPUT_FOLDER/"$name"_temp/1.bam;

PicardCommandLine MarkDuplicates REMOVE_DUPLICATES=true DUPLICATE_SCORING_STRATEGY=TOTAL_MAPPED_REFERENCE_LENGTH \
	M=$INPUT_FOLDER/"$name"_temp/2.txt \
	I=$INPUT_FOLDER/"$name"_temp/1.bam \
	O=$INPUT_FOLDER/"$name"_temp/2.bam;

PicardCommandLine MarkDuplicates REMOVE_DUPLICATES=true DUPLICATE_SCORING_STRATEGY=RANDOM \
	M=$INPUT_FOLDER/"$name"_temp/3.txt \
	I=$INPUT_FOLDER/"$name"_temp/2.bam \
	O=$INPUT_FOLDER/"$name"_temp/3.bam;

samtools markdup -@ 10 -r -S -s -O BAM $INPUT_FOLDER/"$name"_temp/3.bam $INPUT_FOLDER/"$name"_60M_sorted_CLEARED.bam > $INPUT_FOLDER/"$name"_temp/4.txt;

cat $INPUT_FOLDER/"$name"_temp/*.txt >> $INPUT_FOLDER/"$name"_CLEARED_metrics.txt;

rm -rf $INPUT_FOLDER/"$name"_temp;

done
