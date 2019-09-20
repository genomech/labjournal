#!/bin/bash

THREADS=8
B2_INDEX="/dev/datasets/FairWind/_db/hg19/hg19.fa"
OUTPUT_PATH="/dev/datasets/FairWind/_results/bowtie/vcf"
SORTED_PATH="/dev/datasets/FairWind/_results/bowtie/bam"
DEPTH=10
QUALITY=30

for var in '1-1' '1-2' '1-5' '1-6' '1-9'
do
bcftools mpileup --threads $THREADS -f $B2_INDEX $SORTED_PATH/sample-"$var"_sorted.bam | bcftools call --threads $THREADS -cv -Ou | bcftools filter -i "DP>$DEPTH & %QUAL>$QUALITY" > $OUTPUT_PATH/filterCalls_"$var"_DP"$DEPTH"_QUAL"$QUALITY".txt
echo File "$var" is ready.
done
