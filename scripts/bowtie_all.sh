#!/bin/bash

THREADS=12
B2_INDEX="/dev/datasets/FairWind/_db/hg19/hg19_small/hg19"
INPUT_PATH="/dev/datasets/FairWind/_results/cut/splitted"
OUTPUT_PATH="/dev/datasets/FairWind/_results/bowtie/sam"
SORTED_PATH="/dev/datasets/FairWind/_results/bowtie/bam"

mkdir $SORTED_PATH
mkdir $OUTPUT_PATH

for var in '1-1' '1-2' '1-4' '1-5' '1-6' '1-9'
do

bowtie2 -x $B2_INDEX \
 	-U $INPUT_PATH/sample-"$var"_R1_SplitLeft.fastq.gz \
 	-U $INPUT_PATH/sample-"$var"_R1_SplitRight.fastq.gz \
 	-U $INPUT_PATH/sample-"$var"_R2_SplitLeft.fastq.gz \
 	-U $INPUT_PATH/sample-"$var"_R2_SplitRight.fastq.gz \
 	--very-sensitive -p $THREADS | tee $OUTPUT_PATH/sample-"$var".sam | samtools view -bS \
 	-@ $THREADS - | samtools sort -@ $THREADS -O BAM - > $SORTED_PATH/sample-"$var"_sorted.bam;

echo $var is ready.

done

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

 echo 1-3-4 is ready.


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

 echo 1-7-8 is ready.

md5sum $OUTPUT_PATH/*.sam > $OUTPUT_PATH/all.md5;
md5sum $SORTED_PATH/*.bam > $SORTED_PATH/all.md5;
