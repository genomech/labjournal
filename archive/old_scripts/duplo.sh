#!/bin/bash

INPUT_FOLDER="/dev/datasets/FairWind/_results/cut/uncut_bam";
OUTPUT_FOLDER="/dev/datasets/FairWind/_results/cut/uncut_bam_dupless";
mkdir -p $OUTPUT_FOLDER;
for file in $INPUT_FOLDER/*.bam;
do {
start_time=$(StartTime);
PicardCommandLine MarkDuplicates REMOVE_DUPLICATES=true M=$OUTPUT_FOLDER/$(FileBase $file)_metrics.txt I=$file O=$OUTPUT_FOLDER/$(FileBase $file)_dupless.bam;
echo "Sample "$(FileBase $file)" is ready "$(Timestamp $start_time)"";
} done
