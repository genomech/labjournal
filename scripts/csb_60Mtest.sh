#!/bin/bash

THREADS=10

INPUT_FOLDER='/dev/datasets/FairWind/_results/cut/bridgeless_SL'
OUTPUT_FOLDER='/dev/datasets/FairWind/_results/cut/splitted'

mkdir $OUTPUT_FOLDER;

bridge='GCTGAGG'
egdirb='CCTCAGC'
gatc='GATC'

min_length=8
error=0.2

## SPLIT

for var in '1-1' '1-2' '1-3' '1-4' '1-5' '1-6' '1-7' '1-8' '1-9'
do

split_seq=$(if [ "$var" = "1-4" ]; then (echo "print(\""$gatc""$gatc"\")" | python3); else (echo "print(\""$bridge""$gatc""$egdirb"\")" | python3); fi)
overlap=$(echo "print(len(\""$split_seq"\") - 2)" | python3)
echo "*" Start file $var with sequence $split_seq and overlap $overlap...

# 1. Left

cutadapt -m $min_length -e $error -O $overlap -j $THREADS -y ":SPLIT-LEFT:R1" -a "$split_seq" \
    -o $OUTPUT_FOLDER/sample-"$var"_R1_SplitLeft.fastq.gz \
       $INPUT_FOLDER/sample-"$var"_R1_Bridgeless.fastq.gz > $OUTPUT_FOLDER/sample-"$var"_report.txt;
       
cutadapt -m $min_length -e $error -O $overlap -j $THREADS -y ":SPLIT-LEFT:R2" -a "$split_seq" \
    -o $OUTPUT_FOLDER/sample-"$var"_R2_SplitLeft.fastq.gz \
       $INPUT_FOLDER/sample-"$var"_R2_Bridgeless.fastq.gz >> $OUTPUT_FOLDER/sample-"$var"_report.txt;
       
# 2. Right (trimmed-only)

cutadapt -m $min_length -e $error -O $overlap -j $THREADS -y ":SPLIT-RIGHT:R1" -g "$split_seq" --trimmed-only \
    -o $OUTPUT_FOLDER/sample-"$var"_R1_SplitRight.fastq.gz \
       $INPUT_FOLDER/sample-"$var"_R1_Bridgeless.fastq.gz >> $OUTPUT_FOLDER/sample-"$var"_report.txt;
       
cutadapt -m $min_length -e $error -O $overlap -j $THREADS -y ":SPLIT-RIGHT:R2" -g "$split_seq" --trimmed-only \
    -o $OUTPUT_FOLDER/sample-"$var"_R2_SplitRight.fastq.gz \
       $INPUT_FOLDER/sample-"$var"_R2_Bridgeless.fastq.gz >> $OUTPUT_FOLDER/sample-"$var"_report.txt;
       
echo "*" Done.
done
