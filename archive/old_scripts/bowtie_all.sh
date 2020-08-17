#!/bin/bash

function timestamp { PASSED=$(echo $(date +%s) - $2 | bc); echo "* Process '"$1"' done ["$(date -d@$PASSED -u +%H:%M:%S)"]"; }
function seal_folder { md5deep -lr $1/$2 > $1/all.md5; chmod -R 555 $1; echo "["$(date +'%Y-%m-%d %H:%M:%S')"] Folder "$1" is sealed."; }

THREADS=10
B2_INDEX="/dev/datasets/FairWind/_db/hg19/hg19_small/hg19"
INPUT_PATH="/dev/datasets/FairWind/_results/cut/splitted"
OUTPUT_PATH="/dev/datasets/FairWind/_results/cut/sam"
SORTED_PATH="/dev/datasets/FairWind/_results/cut/bam"

mkdir -p $SORTED_PATH
mkdir -p $OUTPUT_PATH

for var in '1-1' '1-2' '1-5' '1-6' '1-9'
do

start_time=$(date +%s);

bowtie2 -x $B2_INDEX \
	-U $INPUT_PATH/sample-"$var"_R1_SplitLeft.fastq.gz \
	-U $INPUT_PATH/sample-"$var"_R1_SplitRight.fastq.gz \
	-U $INPUT_PATH/sample-"$var"_R2_SplitLeft.fastq.gz \
	-U $INPUT_PATH/sample-"$var"_R2_SplitRight.fastq.gz \
	--very-sensitive -p $THREADS | tee $OUTPUT_PATH/sample-"$var".sam | samtools view -bS \
	-@ $THREADS - | samtools sort -@ $THREADS -O BAM - > $SORTED_PATH/sample-"$var"_sorted.bam;

timestamp sample-"$var" $start_time;
	
done

start_time=$(date +%s);

bowtie2 -x $B2_INDEX \
	-U $INPUT_PATH/sample-1-3_R1_SplitLeft.fastq.gz \
	-U $INPUT_PATH/sample-1-3_R1_SplitRight.fastq.gz \
	-U $INPUT_PATH/sample-1-3_R2_SplitLeft.fastq.gz \
	-U $INPUT_PATH/sample-1-3_R2_SplitRight.fastq.gz \
	-U $INPUT_PATH/sample-1-4_R1_SplitLeft.fastq.gz \
	-U $INPUT_PATH/sample-1-4_R1_SplitRight.fastq.gz \
	-U $INPUT_PATH/sample-1-4_R2_SplitLeft.fastq.gz \
	-U $INPUT_PATH/sample-1-4_R2_SplitRight.fastq.gz \
	--very-sensitive -p $THREADS | tee $OUTPUT_PATH/sample-1-3-4.sam | samtools view -bS \
	-@ $THREADS - | samtools sort -@ $THREADS -O BAM - > $SORTED_PATH/sample-1-3-4_sorted.bam;

timestamp sample-1-3-4 $start_time;

start_time=$(date +%s);

bowtie2 -x $B2_INDEX \
	-U $INPUT_PATH/sample-1-7_R1_SplitLeft.fastq.gz \
	-U $INPUT_PATH/sample-1-7_R1_SplitRight.fastq.gz \
	-U $INPUT_PATH/sample-1-7_R2_SplitLeft.fastq.gz \
	-U $INPUT_PATH/sample-1-7_R2_SplitRight.fastq.gz \
	-U $INPUT_PATH/sample-1-8_R1_SplitLeft.fastq.gz \
	-U $INPUT_PATH/sample-1-8_R1_SplitRight.fastq.gz \
	-U $INPUT_PATH/sample-1-8_R2_SplitLeft.fastq.gz \
	-U $INPUT_PATH/sample-1-8_R2_SplitRight.fastq.gz \
	--very-sensitive -p $THREADS | tee $OUTPUT_PATH/sample-1-7-8.sam | samtools view -bS \
	-@ $THREADS - | samtools sort -@ $THREADS -O BAM - > $SORTED_PATH/sample-1-7-8_sorted.bam;

timestamp sample-1-7-8 $start_time;

start_time=$(date +%s);

seal_folder "$OUTPUT_PATH" "*.sam";
seal_folder "$SORTED_PATH" "*.bam";

timestamp SEAL $start_time;
