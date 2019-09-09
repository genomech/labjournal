#!/bin/bash


FAI="/dev/datasets/FairWind/_db/hg19/hg19.fa.fai"
EXOME="/dev/datasets/FairWind/_db/MedExome_hg19_capture_targets.sorted.bed"
INPUT_PATH="/dev/datasets/FairWind/_results/60m/bowtie/bam"
OUTPUT_PATH="/dev/datasets/FairWind/_results/60m/coverage"

mkdir $OUTPUT_PATH

for var in '1-1' '1-2' '1-3' '1-4' '1-5' '1-6' '1-7' '1-8' '1-9'
do

bedtools coverage -hist -sorted \
	-g $FAI \
	-a $EXOME \
	-b $INPUT_PATH/sample-"$var"_sorted.bam > $OUTPUT_PATH/sample-"$var"_ExomeCoverage.txt

echo $var is ready.

done

md5sum $OUTPUT_PATH/*.txt > $OUTPUT_PATH/all.md5;
