#!/bin/bash

# Index Ref

Ref="/dev/datasets/FairWind/_results/NarimanBarcodes/Plasmida/PlasmidaNariman_Barcode.fa"; samtools faidx $Ref; bwa index $Ref; gatk CreateSequenceDictionary -R $Ref;

# Align

Threads=10;
R1="/dev/datasets/ngs_data/Battulin_BGI_20200918/demult/200915_X578_FCHCMFGCCX2_L8_CHKPE85220080212_1_single_ASA23.fq.gz";
R2="/dev/datasets/ngs_data/Battulin_BGI_20200918/demult/200915_X578_FCHCMFGCCX2_L8_CHKPE85220080212_2_single_ASA23.fq.gz";
QueryBam="/dev/datasets/FairWind/_results/NarimanBarcodes/200915_X578_FCHCMFGCCX2_L8_CHKPE85220080212_single_ASA23/200915_X578_FCHCMFGCCX2_L8_CHKPE85220080212_single_ASA23.bam";
StatsTxt="/dev/datasets/FairWind/_results/NarimanBarcodes/200915_X578_FCHCMFGCCX2_L8_CHKPE85220080212_single_ASA23/flagstat.txt";

bwa mem -t $Threads $Ref $R1 $R2 | grep -Pv '[\t]{2,}' | gatk SortSam --VERBOSITY ERROR -SO queryname -I /dev/stdin -O $QueryBam; samtools flagstat $QueryBam > $StatsTxt;

# Slice for IGV

SliceSize=100000;
CoordBam="/dev/datasets/FairWind/_results/NarimanBarcodes/200915_X578_FCHCMFGCCX2_L8_CHKPE85220080212_single_ASA23/slice.bam";

samtools view -h $QueryBam | head -l $SliceSize | gatk SortSam --VERBOSITY ERROR -SO coordinate -I /dev/stdin -O $CoordBam; samtools index $CoordBam;
