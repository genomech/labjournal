#!/bin/bash

QUALIMAP_FOLDER="/dev/datasets/FairWind/_tools/qualimap_v2.2.1"
INPUT_FOLDER="/dev/datasets/FairWind/_results/60m/PRIMARY_ANALYSIS_13D/bam_sorted"
REPORT_FOLDER="/dev/datasets/FairWind/_results/60m/PRIMARY_ANALYSIS_13D/qualimap_reports"
EXOME="/dev/datasets/FairWind/_db/MedExome_hg19_capture_targets_sorted_FOR_HIS_MAGESTY_QUALIMAP.bed"

mkdir -l $REPORT_FOLDER

for name in 'dinara_38_S4_60M_sorted' 'dinara_38_S4_60M_sorted_dupless' 'sample-1-1_60M_sorted' 'sample-1-1_60M_sorted_dupless' 'sample-1-3_60M_sorted' 'sample-1-3_60M_sorted_dupless'
do

# SKIP DUPLICATES
$QUALIMAP_FOLDER/qualimap bamqc -bam $INPUT_FOLDER/"$name".bam --java-mem-size=20G \
	-c --feature-file $EXOME --outside-stats \
	-outdir $REPORT_FOLDER -outfile "$name"_report_skipdup.pdf -outformat PDF \
	--skip-duplicated --skip-dup-mode 2 ;

# NOT SKIP
$QUALIMAP_FOLDER/qualimap bamqc -bam $INPUT_FOLDER/"$name".bam --java-mem-size=20G \
	-c --feature-file $EXOME --outside-stats \
	-outdir $REPORT_FOLDER -outfile "$name"_report_notskip.pdf -outformat PDF ;

echo "$name" is ready.

done 
