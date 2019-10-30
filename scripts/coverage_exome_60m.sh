#!/bin/bash

FAI="/dev/datasets/FairWind/_db/hg19/hg19.fa.fai"
EXOME="/dev/datasets/FairWind/_db/MedExome_hg19_capture_targets.sorted.bed"
INPUT_PATH="/dev/datasets/FairWind/_results/"
OUTPUT_PATH="/dev/datasets/FairWind/_results/"

mkdir $OUTPUT_PATH

for var in '38_S4'
do

bedtools coverage -hist -sorted \
	-g $FAI \
	-a $EXOME \
	-b $INPUT_PATH/"$var".bam > $INPUT_PATH/"$var"--CoverageExome.txt

echo $var is ready.

done

