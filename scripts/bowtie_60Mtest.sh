#!/bin/bash

THREADS=10
B2_INDEX="/dev/datasets/FairWind/_db/hg19/hg19_small/hg19"
INPUT_PATH="/dev/datasets/FairWind/_results/60m/PRIMARY_ANALYSIS_13D/ready_fastq"
OUTPUT_PATH="/dev/datasets/FairWind/_results/60m/PRIMARY_ANALYSIS_13D/sam"
SORTED_PATH="/dev/datasets/FairWind/_results/60m/PRIMARY_ANALYSIS_13D/bam_sorted"

mkdir $SORTED_PATH
mkdir $OUTPUT_PATH

for var in '1-1' '1-3'
do
bowtie2 -x $B2_INDEX \
 	-U $INPUT_PATH/sample-"$var"_R1_SplitLeft.fastq.gz \
 	-U $INPUT_PATH/sample-"$var"_R1_SplitRight.fastq.gz \
 	-U $INPUT_PATH/sample-"$var"_R2_SplitLeft.fastq.gz \
 	-U $INPUT_PATH/sample-"$var"_R2_SplitRight.fastq.gz \
 	--very-sensitive -p $THREADS | tee $OUTPUT_PATH/sample-"$var"_60M.sam | samtools view -bS \
 	-@ $THREADS - | samtools sort -@ $THREADS -O BAM - > $SORTED_PATH/sample-"$var"_60M_sorted.bam;

echo $var is ready.
done

# bowtie2 -x $B2_INDEX \
# 	-1 $INPUT_PATH/dinara_38_S4_R1_60M_Illuminaless.fastq.gz \
# 	-2 $INPUT_PATH/dinara_38_S4_R2_60M_Illuminaless.fastq.gz \
# 	--very-sensitive -p $THREADS | tee $OUTPUT_PATH/dinara_38_S4_60M.sam | samtools view -bS \
# 	-@ $THREADS - | samtools sort -@ $THREADS -O BAM - > $SORTED_PATH/dinara_38_S4_60M_sorted.bam;
# 
# echo Dinara is ready.

md5sum $OUTPUT_PATH/*.sam > $OUTPUT_PATH/all.md5;
md5sum $SORTED_PATH/*.bam > $SORTED_PATH/all.md5;
