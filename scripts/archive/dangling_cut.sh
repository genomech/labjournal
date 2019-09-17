 #!/bin/bash

THREADS=12

INPUT_FOLDER='/dev/datasets/FairWind/_results/dangling_ends'
OUTPUT_FOLDER='/dev/datasets/FairWind/_results/dangling_cut'

mkdir $OUTPUT_FOLDER;

## CUT ILLUMINA

illumina1='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
illumina2='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'

for var in '3' '5' '7' '8'
do
(
cutadapt -m 8 -j $THREADS -a $illumina1 -A $illumina2 \
    -o $OUTPUT_FOLDER/s"$var"_FR_RF_0_R1_Illuminaless.fastq.gz \
    -p $OUTPUT_FOLDER/s"$var"_FR_RF_0_R2_Illuminaless.fastq.gz \
       $INPUT_FOLDER/s"$var"_FR_RF_0_R1.fastq.gz \
       $INPUT_FOLDER/s"$var"_FR_RF_0_R2.fastq.gz; ) > $OUTPUT_FOLDER/s"$var"_FR_RF_0_cutadapt.txt;

done

md5sum $OUTPUT_FOLDER/*.fastq.gz > $OUTPUT_FOLDER/illuminaless/all.md5;
