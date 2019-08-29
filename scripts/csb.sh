#!/bin/bash

THREADS=12

INPUT_FOLDER='/dev/datasets/ngs_data/ExoC_Belopuz/30-213832944'
OUTPUT_FOLDER='/dev/datasets/FairWind/_results/cut_and_split'

B2_INDEX="/dev/datasets/FairWind/_db/hg19/hg19_small/hg19"

mkdir $OUTPUT_FOLDER;
mkdir $OUTPUT_FOLDER/cut_reports;
mkdir $OUTPUT_FOLDER/splitted;
mkdir $OUTPUT_FOLDER/sam;
mkdir $OUTPUT_FOLDER/bam;

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

for var in '1-4'
do
(echo "---------------- THIS IS $var (SPARTAAAAA) !!!! ----------------";

## CUT ILLUMINA

cutadapt -m 8 -j $THREADS -a $illumina1 -A $illumina2 \
    -o $OUTPUT_FOLDER/sample-"$var"_R1_TEMP0.fastq \
    -p $OUTPUT_FOLDER/sample-"$var"_R2_TEMP0.fastq \
       $INPUT_FOLDER/sample-"$var"_R1_001.fastq.gz \
       $INPUT_FOLDER/sample-"$var"_R2_001.fastq.gz;

## CUT BRIDGE ENDS

cutadapt -m 8 -e 0.2 -O 7 -j $THREADS -a "$double_seq"X -g X"$double_seq" -A "$double_seq"X -G X"$double_seq" \
    -o $OUTPUT_FOLDER/sample-"$var"_R1_TEMP1.fastq \
    -p $OUTPUT_FOLDER/sample-"$var"_R2_TEMP1.fastq \
       $OUTPUT_FOLDER/sample-"$var"_R1_TEMP0.fastq \
       $OUTPUT_FOLDER/sample-"$var"_R2_TEMP0.fastq;

cutadapt -m 8 -e 0.2 -O 7 -j $THREADS -a "$double_2_right"X -g X"$double_2_left" -A "$double_2_right"X -G X"$double_2_left" \
    -o $OUTPUT_FOLDER/sample-"$var"_R1_TEMP0.fastq \
    -p $OUTPUT_FOLDER/sample-"$var"_R2_TEMP0.fastq \
       $OUTPUT_FOLDER/sample-"$var"_R1_TEMP1.fastq \
       $OUTPUT_FOLDER/sample-"$var"_R2_TEMP1.fastq;

## SPLIT

# 1. Left

cutadapt -m 8 -e 0.2 -O 7 -j $THREADS -y SPLIT_LEFT -a "$split_seq" \
    -o $OUTPUT_FOLDER/splitted/sample-"$var"_R1_SplitLeft.fastq \
       $OUTPUT_FOLDER/sample-"$var"_R1_TEMP0.fastq;
       
cutadapt -m 8 -e 0.2 -O 7 -j $THREADS -y SPLIT_LEFT -a "$split_seq" \
    -o $OUTPUT_FOLDER/splitted/sample-"$var"_R2_SplitLeft.fastq \
       $OUTPUT_FOLDER/sample-"$var"_R2_TEMP0.fastq;
       
# 2. Right (trimmed-only)

cutadapt -m 8 -e 0.2 -O 7 --trimmed-only -y SPLIT_RIGHT -j $THREADS -g "$split_seq" \
    -o $OUTPUT_FOLDER/splitted/sample-"$var"_R1_SplitRight.fastq \
       $OUTPUT_FOLDER/sample-"$var"_R1_TEMP0.fastq;
       
cutadapt -m 8 -e 0.2 -O 7 --trimmed-only -y SPLIT_RIGHT -j $THREADS -g "$split_seq" \
    -o $OUTPUT_FOLDER/splitted/sample-"$var"_R2_SplitRight.fastq \
       $OUTPUT_FOLDER/sample-"$var"_R2_TEMP0.fastq;

# 3. Remove temp files

rm -rf -v $OUTPUT_FOLDER/sample-"$var"_R*_TEMP*.fastq;

# 4. Bowtie

#bowtie2 -x $B2_INDEX \
#	-U $OUTPUT_FOLDER/splitted/sample-"$var"_R1_SplitLeft.fastq \
#	-U $OUTPUT_FOLDER/splitted/sample-"$var"_R1_SplitRight.fastq \
#	-U $OUTPUT_FOLDER/splitted/sample-"$var"_R2_SplitLeft.fastq \
#	-U $OUTPUT_FOLDER/splitted/sample-"$var"_R2_SplitRight.fastq \
#	--very-sensitive -p $THREADS | tee $OUTPUT_FOLDER/sam/sample-"$var".sam | samtools view -bS \
#	-@ $THREADS - | samtools sort -@ $THREADS -O BAM - > $OUTPUT_FOLDER/bam/sample-"$var"_sorted.bam \
) > $OUTPUT_FOLDER/cut_reports/sample-"$var"_report.txt;

done
