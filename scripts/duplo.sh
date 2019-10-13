#!/bin/bash

INPUT_FOLDER="/dev/datasets/FairWind/_results/bowtie/bam"
OUTPUT_FOLDER="/dev/datasets/FairWind/_results/bowtie/dupless"

mkdir -l $OUTPUT_FOLDER

for name in '1-1' '1-2' '1-3-4' '1-5' '1-6' '1-7-8' '1-9'
do
PicardCommandLine MarkDuplicates REMOVE_DUPLICATES=true \
	M=$OUTPUT_FOLDER/sample-"$name"_metrics.txt \
	I=$INPUT_FOLDER/sample-"$name"_sorted.bam \
	O=$OUTPUT_FOLDER/sample-"$name"_dupless.bam;
done
