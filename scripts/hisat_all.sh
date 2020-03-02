#!/bin/bash

THREADS=12
HISAT2_INDEX="/dev/datasets/FairWind/_db/rn6/rn6"
INPUT_PATH="/dev/datasets/ngs_data/SeqDataTransfer/AndreFelipeRodrigues"
SAM_PATH="/dev/datasets/FairWind/_results/Andre/sam"
BAM_PATH="/dev/datasets/FairWind/_results/Andre/bam_sorted"

mkdir -p $SAM_PATH
mkdir -p $BAM_PATH

for var in '001_S1' '002_S2' '003_S3' '004_S4' '005_S5' '006_S6' '007_S7' '008_S8' '009_S9' '010_S10' '011_S11' '012_S12'
do
hisat2 -x $HISAT2_INDEX \
	-1 $INPUT_PATH/MB_AFR_"$var"_R1_001.fastq.gz \
	-2 $INPUT_PATH/MB_AFR_"$var"_R1_001.fastq.gz \
	-p $THREADS | tee $SAM_PATH/MB_AFR_"$var".sam | samtools view -bS \
	-@ $THREADS - | samtools sort -@ $THREADS -O BAM - > $BAM_PATH/MB_AFR_"$var"_sorted.bam;
echo $var is ready.
done

md5sum $SAM_PATH/*.sam > $SAM_PATH/all.md5;
chmod 555 -R $SAM_PATH;
md5sum $BAM_PATH/*.bam > $BAM_PATH/all.md5;
chmod 555 -R $BAM_PATH;

echo Dirs are sealed.
