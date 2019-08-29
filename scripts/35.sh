#!/bin/bash

THREADS=12

INPUT_FOLDER='/dev/datasets/ngs_data/ExoC_Belopuz/30-213832944'
OUTPUT_FOLDER='/dev/datasets/FairWind/_results/35'

B2_INDEX="/dev/datasets/FairWind/_db/hg19/hg19_small/hg19"

mkdir $OUTPUT_FOLDER;
mkdir $OUTPUT_FOLDER/cut_reports;

illumina1='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
illumina2='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'

bridge='GCTGAGG'
egdirb='CCTCAGC'
gatc='GATC'

double_seq="$bridge""$gatc""$egdirb""$bridge""$gatc""$egdirb"

double_2_left="$bridge""$gatc""$egdirb""$bridge"
double_2_right="$egdirb""$bridge""$gatc""$egdirb"

split_seq="$gatc""$gatc"
#split_seq="$bridge""$gatc""$egdirb"

for var in '1-3' '1-5'
do
(echo "---------------- THIS IS $var (SPARTAAAAA) !!!! ----------------";

## CUT ILLUMINA

cutadapt -m 8 -j $THREADS -a $illumina1 -A $illumina2 \
    -o $OUTPUT_FOLDER/sample-"$var"_R1_Illuminaless.fastq \
    -p $OUTPUT_FOLDER/sample-"$var"_R2_Illuminaless.fastq \
       $INPUT_FOLDER/sample-"$var"_R1_001.fastq.gz \
       $INPUT_FOLDER/sample-"$var"_R2_001.fastq.gz;
) > $OUTPUT_FOLDER/cut_reports/sample-"$var"_report.txt;

done
