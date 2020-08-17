#!/bin/bash

THREADS=12
B2_INDEX="/dev/datasets/FairWind/hg19_small/hg19"
INPUT_PATH="/dev/datasets/FairWind/cut_stuff"
OUTPUT_PATH="./sams"
SORTED_PATH="./bams"
METRICS_PATH="./metrics"
ERR_PATH="./err"

# 1-4 alone

( bowtie2 -x $B2_INDEX \
    -U $INPUT_PATH/sample-1-4_R1_SplitLeft.fastq.gz \
    -U $INPUT_PATH/sample-1-4_R2_SplitLeft.fastq.gz \
    -U $INPUT_PATH/sample-1-4_R1_SplitRight.fastq.gz \
    -U $INPUT_PATH/sample-1-4_R2_SplitRight.fastq.gz \
    --very-sensitive -p $THREADS --met-file $METRICS_PATH/sample-1-4_metrics.txt | tee $OUTPUT_PATH/sample-1-4.sam | samtools view -bS -@ $THREADS - | samtools sort -@ $THREADS -O BAM - > $SORTED_PATH/sample-1-4.sorted.bam ) 3>&1 1>&2 2>&3 | tee $ERR_PATH/sample-1-4_stderr.log;
    
echo File sample-1-4 is ready.
