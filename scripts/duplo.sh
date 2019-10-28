#!/bin/bash

INPUT_FOLDER="/dev/datasets/FairWind/_results/60m/PRIMARY_ANALYSIS_13D/bam_sorted"

for name in 'sample-1-1'
do

samtools markdup -@ 10 -r -S -s -O BAM $INPUT_FOLDER/"$name"_60M_sorted.bam $INPUT_FOLDER/"$name"_60M_sorted_dupless_SAMTOOLS.bam

#PicardCommandLine MarkDuplicates REMOVE_DUPLICATES=true DUPLICATE_SCORING_STRATEGY=TOTAL_MAPPED_REFERENCE_LENGTH \
#	M=$INPUT_FOLDER/"$name"_dupless_metrics_ANOTHER.txt \
#	I=$INPUT_FOLDER/"$name"_60M_sorted.bam \
#	O=$INPUT_FOLDER/"$name"_60M_sorted_dupless_ANOTHER.bam;

done
