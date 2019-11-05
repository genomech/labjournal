#!/bin/bash

FAI="/dev/datasets/FairWind/_db/hg19/hg19.fa.fai"
GENOME="/dev/datasets/FairWind/_db/hg19/hg19.bed"
NOT_EXOME="/dev/datasets/FairWind/_db/NOT_MedExome_hg19.bed"
INPUT_PATH="/dev/datasets/FairWind/_results/60m/PRIMARY_ANALYSIS_13D/bam_sorted"
OUTPUT_PATH="/dev/datasets/FairWind/_results/60m/PRIMARY_ANALYSIS_13D/coverage/full"

mkdir -l $OUTPUT_PATH;

for var in 'dinara_38_S4_60M_sorted' 'dinara_38_S4_60M_sorted_dupless' 'sample-1-1_60M_sorted' 'sample-1-1_60M_sorted_dupless_2' 'sample-1-1_60M_sorted_dupless' 'sample-1-3_60M_sorted' 'sample-1-3_60M_sorted_dupless' 
do

bedtools coverage -hist -sorted \
	-g $FAI \
    -a $GENOME \
    -b $INPUT_PATH/"$var".bam > $OUTPUT_PATH/"$var"_GenomeCoverage.txt;
    
echo $var Genome is ready.

bedtools coverage -hist -sorted \
	-g $FAI \
	-a $NOT_EXOME \
	-b $INPUT_PATH/"$var".bam > $OUTPUT_PATH/"$var"_NotExomeCoverage.txt;

echo $var NotExome is ready.

done
