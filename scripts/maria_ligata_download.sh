#!/bin/bash

OUTPUT_PATH="/dev/datasets/ngs_data/ma_add"

for code in SRR5063167 SRR2033054
do
fastq-dump -Z --split-3 $code | head -n 4000000 | gzip -c - > $OUTPUT_PATH/"$code"_1M.fastq.gz;
echo $code is done.
done

