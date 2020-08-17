#!/bin/bash

INPUT_FOLDER='/dev/datasets/ngs_data/biblexome'
OUTPUT_FOLDER='/dev/datasets/FairWind/_results/60m_dinara'

mkdir $OUTPUT_FOLDER;

bridge='GCTGAGG'
egdirb='CCTCAGC'
gatc='GATC'

## CUT 60M

mkdir $OUTPUT_FOLDER/samples;

for var in '104_S3' '111_S6' '113_S2' '117_S5' '38_S4' '98_S1' 'le1_S7' 'le2_S8' 'le3_S9' 'le4_S10' 'le5_S11' 'le6_S12'
do

gzip -cd $INPUT_FOLDER/"$var"_R1_001.fastq.gz | head -240000000 > $OUTPUT_FOLDER/samples/"$var"_R1_60M.fastq;
gzip -cd $INPUT_FOLDER/"$var"_R2_001.fastq.gz | head -240000000 > $OUTPUT_FOLDER/samples/"$var"_R2_60M.fastq;

done

md5sum $OUTPUT_FOLDER/samples/*.fastq > $OUTPUT_FOLDER/samples/all.md5;

echo Samples are done.
