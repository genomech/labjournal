#!/bin/bash

THREADS=12

INPUT_FOLDER='/dev/datasets/ngs_data/ExoC_Belopuz/30-213832944'
OUTPUT_FOLDER='/dev/datasets/FairWind/_results/cut'

mkdir $OUTPUT_FOLDER;

bridge='GCTGAGG'
egdirb='CCTCAGC'
gatc='GATC'

## CUT ILLUMINA

illumina1='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
illumina2='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'

mkdir $OUTPUT_FOLDER/illuminaless;

for var in '1-1' '1-2' '1-3' '1-4' '1-5' '1-6' '1-7' '1-8' '1-9'
do
(
cutadapt -m 8 -j $THREADS -a $illumina1 -A $illumina2 \
    -o $OUTPUT_FOLDER/illuminaless/sample-"$var"_R1_Illuminaless.fastq.gz \
    -p $OUTPUT_FOLDER/illuminaless/sample-"$var"_R2_Illuminaless.fastq.gz \
       $INPUT_FOLDER/sample-"$var"_R1_001.fastq.gz \
       $INPUT_FOLDER/sample-"$var"_R2_001.fastq.gz; ) > $OUTPUT_FOLDER/illuminaless/sample-"$var"_cutadapt.txt;

done

md5sum $OUTPUT_FOLDER/illuminaless/*.fastq.gz > $OUTPUT_FOLDER/illuminaless/all.md5;

echo Illumina is done.

## CUT BRIDGE ENDS

double_seq="$bridge""$gatc""$egdirb""$bridge""$gatc""$egdirb"
double_2_left="$bridge""$gatc""$egdirb""$bridge"
double_2_right="$egdirb""$bridge""$gatc""$egdirb"

mkdir $OUTPUT_FOLDER/bridgeless;

for var in '1-1' '1-2' '1-3' '1-4' '1-5' '1-6' '1-7' '1-8' '1-9'
do
(
cutadapt -m 8 -e 0.2 -O 7 -j $THREADS -a "$double_seq"X -g X"$double_seq" -A "$double_seq"X -G X"$double_seq" \
    -o $OUTPUT_FOLDER/sample-"$var"_R1_TEMP1.fastq \
    -p $OUTPUT_FOLDER/sample-"$var"_R2_TEMP1.fastq \
       $OUTPUT_FOLDER/illuminaless/sample-"$var"_R1_Illuminaless.fastq.gz \
       $OUTPUT_FOLDER/illuminaless/sample-"$var"_R2_Illuminaless.fastq.gz;
 
cutadapt -m 8 -e 0.2 -O 7 -j $THREADS -a "$double_2_right"X -g X"$double_2_left" -A "$double_2_right"X -G X"$double_2_left" \
    -o $OUTPUT_FOLDER/bridgeless/sample-"$var"_R1_Bridgeless.fastq.gz \
    -p $OUTPUT_FOLDER/bridgeless/sample-"$var"_R2_Bridgeless.fastq.gz \
       $OUTPUT_FOLDER/sample-"$var"_R1_TEMP1.fastq \
       $OUTPUT_FOLDER/sample-"$var"_R2_TEMP1.fastq;
       
rm -rf -v $OUTPUT_FOLDER/sample-"$var"_R*_TEMP*.fastq; ) > $OUTPUT_FOLDER/bridgeless/sample-"$var"_cutadapt.txt;
        
done

md5sum $OUTPUT_FOLDER/bridgeless/*.fastq.gz > $OUTPUT_FOLDER/bridgeless/all.md5;

echo Bridge is done.
 
## SPLIT

split_seq4="$gatc""$gatc"
split_seq="$bridge""$gatc""$egdirb"

mkdir $OUTPUT_FOLDER/splitted;

for var in '1-1' '1-2' '1-3' '1-5' '1-6' '1-7' '1-8' '1-9'
do
(
# 1. Left

cutadapt -m 8 -e 0.2 -O 7 -j $THREADS -y SPLIT_LEFT -a "$split_seq" \
    -o $OUTPUT_FOLDER/splitted/sample-"$var"_R1_SplitLeft.fastq.gz \
       $OUTPUT_FOLDER/bridgeless/sample-"$var"_R1_Bridgeless.fastq.gz;
       
cutadapt -m 8 -e 0.2 -O 7 -j $THREADS -y SPLIT_LEFT -a "$split_seq" \
    -o $OUTPUT_FOLDER/splitted/sample-"$var"_R2_SplitLeft.fastq.gz \
       $OUTPUT_FOLDER/bridgeless/sample-"$var"_R2_Bridgeless.fastq.gz;
       
# 2. Right (trimmed-only)

cutadapt -m 8 -e 0.2 -O 7 --trimmed-only -y SPLIT_RIGHT -j $THREADS -g "$split_seq" \
    -o $OUTPUT_FOLDER/splitted/sample-"$var"_R1_SplitRight.fastq.gz \
       $OUTPUT_FOLDER/bridgeless/sample-"$var"_R1_Bridgeless.fastq.gz;
       
cutadapt -m 8 -e 0.2 -O 7 --trimmed-only -y SPLIT_RIGHT -j $THREADS -g "$split_seq" \
    -o $OUTPUT_FOLDER/splitted/sample-"$var"_R2_SplitRight.fastq.gz \
       $OUTPUT_FOLDER/bridgeless/sample-"$var"_R2_Bridgeless.fastq.gz; ) > $OUTPUT_FOLDER/splitted/sample-"$var"_cutadapt.txt;

done

# 1. Left 1-4
(
cutadapt -m 8 -e 0.2 -O 7 -j $THREADS -y SPLIT_LEFT -a "$split_seq4" \
    -o $OUTPUT_FOLDER/splitted/sample-1-4_R1_SplitLeft.fastq.gz \
       $OUTPUT_FOLDER/bridgeless/sample-1-4_R1_Bridgeless.fastq.gz;
       
cutadapt -m 8 -e 0.2 -O 7 -j $THREADS -y SPLIT_LEFT -a "$split_seq4" \
    -o $OUTPUT_FOLDER/splitted/sample-1-4_R2_SplitLeft.fastq.gz \
       $OUTPUT_FOLDER/bridgeless/sample-1-4_R2_Bridgeless.fastq.gz;
       
# 2. Right 1-4 (trimmed-only)

cutadapt -m 8 -e 0.2 -O 7 --trimmed-only -y SPLIT_RIGHT -j $THREADS -g "$split_seq4" \
    -o $OUTPUT_FOLDER/splitted/sample-1-4_R1_SplitRight.fastq.gz \
       $OUTPUT_FOLDER/bridgeless/sample-1-4_R1_Bridgeless.fastq.gz;
       
cutadapt -m 8 -e 0.2 -O 7 --trimmed-only -y SPLIT_RIGHT -j $THREADS -g "$split_seq4" \
    -o $OUTPUT_FOLDER/splitted/sample-1-4_R2_SplitRight.fastq.gz \
       $OUTPUT_FOLDER/bridgeless/sample-1-4_R2_Bridgeless.fastq.gz; ) > $OUTPUT_FOLDER/bridgeless/sample-1-4_cutadapt.txt;

done

md5sum $OUTPUT_FOLDER/splitted/*.fastq.gz > $OUTPUT_FOLDER/splitted/all.md5;

echo Split is done.
