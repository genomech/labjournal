#!/bin/bash

THREADS=10

INPUT_FOLDER='/dev/datasets/FairWind/_results/60m/PRIMARY_ANALYSIS_13D/native_fastq'
OUTPUT_FOLDER='/dev/datasets/FairWind/_results/60m/PRIMARY_ANALYSIS_13D/ready_fastq'

mkdir $OUTPUT_FOLDER;

bridge='GCTGAGG'
egdirb='CCTCAGC'
gatc='GATC'

## SPLIT

split_seq4="$gatc""$gatc"
split_seq="$bridge""$gatc""$egdirb"

for var in '1-1' '1-3'
do

# 1. Left

cutadapt -m 8 -e 0.2 -O 7 -j $THREADS -x "SPLIT-LEFT:R1_" -a "$split_seq" \
    -o $OUTPUT_FOLDER/sample-"$var"_R1_SplitLeft.fastq.gz \
       $INPUT_FOLDER/sample-"$var"_R1_Bridgeless.fastq.gz;
       
cutadapt -m 8 -e 0.2 -O 7 -j $THREADS -x "SPLIT-LEFT:R2_" -a "$split_seq" \
    -o $OUTPUT_FOLDER/sample-"$var"_R2_SplitLeft.fastq.gz \
       $INPUT_FOLDER/sample-"$var"_R2_Bridgeless.fastq.gz;
       
# 2. Right (trimmed-only)

cutadapt -m 8 -e 0.2 -O 7 --trimmed-only -x "SPLIT-RIGHT:R1_" -j $THREADS -g "$split_seq" \
    -o $OUTPUT_FOLDER/sample-"$var"_R1_SplitRight.fastq.gz \
       $INPUT_FOLDER/sample-"$var"_R1_Bridgeless.fastq.gz;
       
cutadapt -m 8 -e 0.2 -O 7 --trimmed-only -x "SPLIT-RIGHT:R2_" -j $THREADS -g "$split_seq" \
    -o $OUTPUT_FOLDER/sample-"$var"_R2_SplitRight.fastq.gz \
       $INPUT_FOLDER/sample-"$var"_R2_Bridgeless.fastq.gz;

done

echo Split is done.
