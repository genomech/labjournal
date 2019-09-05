#!/bin/bash

INPUT_FOLDER='/dev/datasets/ngs_data/ExoC_Belopuz/30-213832944'
OUTPUT_FOLDER='/dev/datasets/FairWind/_results/60m'

mkdir $OUTPUT_FOLDER;

bridge='GCTGAGG'
egdirb='CCTCAGC'
gatc='GATC'

## CUT 60M

mkdir $OUTPUT_FOLDER/samples;

for var in '1-1' '1-2' '1-3' '1-4' '1-5' '1-6' '1-7' '1-8' '1-9'
do

gzip -cd $INPUT_FOLDER/sample-"$var"_R1_001.fastq.gz | head -240000000 > $OUTPUT_FOLDER/samples/sample-"$var"_R1_60M.fastq;
gzip -cd $INPUT_FOLDER/sample-"$var"_R2_001.fastq.gz | head -240000000 > $OUTPUT_FOLDER/samples/sample-"$var"_R2_60M.fastq;

done

md5sum $OUTPUT_FOLDER/samples/*.fastq > $OUTPUT_FOLDER/illuminaless/all.md5;

echo Samples are done.
