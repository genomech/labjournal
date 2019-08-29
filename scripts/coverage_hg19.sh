#!/bin/bash

bedtools coverage -d -sorted \
    -g /dev/datasets/FairWind/_db/hg19/hg19.fa.fai \
    -a /dev/datasets/FairWind/_db/hg19/hg19.bed \
    -b /dev/datasets/FairWind/_results/bowtie/bam/sample-1-4_sorted.bam > /dev/datasets/FairWind/_results/bowtie/coverage/1-4_hg19_coverage_each.txt;
echo Ready.
