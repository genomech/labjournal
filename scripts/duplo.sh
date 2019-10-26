#!/bin/bash

INPUT_FOLDER="/dev/datasets/FairWind/_results/60m/PRIMARY_ANALYSIS_13D/bam_sorted"

for name in 'sample-1-1' 'sample-1-3'
do

PicardCommandLine MarkDuplicates REMOVE_DUPLICATES=true \
	M=$INPUT_FOLDER/"$name"_dupless_metrics.txt \
	I=$INPUT_FOLDER/"$name"_60M_sorted.bam \
	O=$INPUT_FOLDER/"$name"_60M_sorted_dupless.bam;

done
