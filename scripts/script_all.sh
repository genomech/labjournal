#!/bin/bash

THREADS=10
B2_INDEX="/dev/datasets/FairWind/hg19_small/hg19"
INPUT_PATH="/dev/datasets/ngs_data/biblexome"
OUTPUT_PATH="/dev/datasets/FairWind/dinara/samwise"
SORTED_PATH="/dev/datasets/FairWind/dinara/sorted"
METRICS_PATH="/dev/datasets/FairWind/dinara/metrics"
ERR_PATH="/dev/datasets/FairWind/dinara/err"

for var in le6_S12
do
( bowtie2 -x $B2_INDEX -1 $INPUT_PATH/"$var"_R1_001.fastq.gz -2 $INPUT_PATH/"$var"_R2_001.fastq.gz \
	--very-sensitive -p $THREADS --met-file $METRICS_PATH/"$var"_metrics.txt | tee $OUTPUT_PATH/"$var".sam | samtools view -bS \
	-@ $THREADS - | samtools sort -@ $THREADS -O BAM - > $SORTED_PATH/"$var".sorted.bam ) 3>&1 1>&2 2>&3 | tee $ERR_PATH/"$var"_stderr.log
echo Files "$var" are ready.
done
