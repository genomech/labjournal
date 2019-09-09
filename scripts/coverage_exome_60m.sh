#!/bin/bash


FAI="/dev/datasets/FairWind/_db/hg19/hg19.fa.fai"
EXOME="/dev/datasets/FairWind/_db/MedExome_hg19_capture_targets.sorted.bed"
INPUT_PATH="/dev/datasets/FairWind/_results/20-120M/sorted"
OUTPUT_PATH="/dev/datasets/FairWind/_results/20-120M/coverage"

mkdir $OUTPUT_PATH

for var in '20' '40' '60' '80' '100' '120'
do

bedtools coverage -hist -sorted \
	-g $FAI \
	-a $EXOME \
	-b $INPUT_PATH/sample-1-2_"$var"_sorted.bam > $OUTPUT_PATH/sample-1-2_"$var".txt

echo $var is ready.

done

md5sum $OUTPUT_PATH/*.txt > $OUTPUT_PATH/all.md5;
