#!/bin/bash

THREADS=6
INPUT_PATH="/dev/datasets/FairWind/_results/20-120M/bam"
OUTPUT_PATH="/dev/datasets/FairWind/_results/20-120M/sorted"

mkdir $OUTPUT_PATH

for var in '20' '40' '60' '80' '100' '120'
do

samtools view -bS -@ $THREADS $INPUT_PATH/sample-1-2_"$var".bam | samtools sort -@ $THREADS -O BAM - > $OUTPUT_PATH/sample-1-2_"$var"_sorted.bam

echo $var is ready.

done
