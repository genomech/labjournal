# Декабрьские данные (лаборатория Воеводы)

Первые 4 баркода.

1. Качество fastq
2. Выравнять (удалять дубликаты, если низкое покрытие)
3. Сравнить мутации в геноме и экзоме

## FastQC

* [Экзом](http://htmlpreview.github.io/?https://github.com/regnveig/labjournal/blob/master/labjournal/FastQC_results/December_Dinara/v300021589_L03_41_1_fastqc.html)
* Геном:
[1](http://htmlpreview.github.io/?https://github.com/regnveig/labjournal/blob/master/labjournal/FastQC_results/December_Dinara/V300021405_L02_1_1_fastqc.html),
[2](http://htmlpreview.github.io/?https://github.com/regnveig/labjournal/blob/master/labjournal/FastQC_results/December_Dinara/V300021405_L02_2_1_fastqc.html),
[3](http://htmlpreview.github.io/?https://github.com/regnveig/labjournal/blob/master/labjournal/FastQC_results/December_Dinara/V300021405_L02_3_1_fastqc.html),
[4](http://htmlpreview.github.io/?https://github.com/regnveig/labjournal/blob/master/labjournal/FastQC_results/December_Dinara/V300021405_L02_4_1_fastqc.html).

## Обрезка адаптеров

```bash
boomer;
Logo "Illuminaless (cutadapt)";
THREADS=10;
OUTPUT_FOLDER='/dev/datasets/FairWind/_results/December_Dinara/Adaptless';
mkdir -p $OUTPUT_FOLDER;
illumina1='AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA';
illumina2='AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG';
r1=(/dev/datasets/ngs_data/Dinara_December/Китай/MGI_helicon/*/*1.fq.gz);
r2=(/dev/datasets/ngs_data/Dinara_December/Китай/MGI_helicon/*/*2.fq.gz);
for var in ${!r1[*]};
do {
start_time=$(StartTime);
r1_output=$OUTPUT_FOLDER/$(FileBase ${r1[var]})_Adaptless.fq.gz;
r2_output=$OUTPUT_FOLDER/$(FileBase ${r2[var]})_Adaptless.fq.gz;
fb=$(FileBase ${r1[var]});
report=$OUTPUT_FOLDER/${fb::-2}_report.txt;
cutadapt -m 8 -j $THREADS -a $illumina1 -A $illumina2 -o $r1_output -p $r2_output ${r1[var]} ${r2[var]} > $report;
echo "File "${fb::-2}" is done "$(Timestamp $start_time)"";
} done;
Seal $OUTPUT_FOLDER
```

И ещё мусор:

```bash
boomer;
Logo "Trashless (cutadapt)";
THREADS=10;
OUTPUT_FOLDER='/dev/datasets/FairWind/_results/December_Dinara/Trashless';
mkdir -p $OUTPUT_FOLDER;
illumina1='TTGTCTTCCTAAGGAACGACATGGCTACGATCCGACTT';
r1=(/dev/datasets/FairWind/_results/December_Dinara/Adaptless/*1_Adaptless.fq.gz);
r2=(/dev/datasets/FairWind/_results/December_Dinara/Adaptless/*2_Adaptless.fq.gz);
for var in ${!r1[*]};
do {
start_time=$(StartTime);
r1_output=$OUTPUT_FOLDER/$(FileBase ${r1[var]})_Trashless.fq.gz;
r2_output=$OUTPUT_FOLDER/$(FileBase ${r2[var]})_Trashless.fq.gz;
fb=$(FileBase ${r1[var]});
report=$OUTPUT_FOLDER/${fb::-2}_report.txt;
cutadapt -m 8 -j $THREADS -a $illumina1 -o $r1_output -p $r2_output ${r1[var]} ${r2[var]} > $report;
echo "File "${fb::-2}" is done "$(Timestamp $start_time)"";
} done;
Seal $OUTPUT_FOLDER
```

Мусора почти не осталось.

* [Экзом](http://htmlpreview.github.io/?https://github.com/regnveig/labjournal/blob/master/labjournal/FastQC_results/Trashless_fqc_December/v300021589_L03_41_1_Adaptless_Trashless_fastqc.html)
* Геном:
[1](http://htmlpreview.github.io/?https://github.com/regnveig/labjournal/blob/master/labjournal/FastQC_results/Trashless_fqc_December/V300021405_L02_1_1_Adaptless_Trashless_fastqc.html),
[2](http://htmlpreview.github.io/?https://github.com/regnveig/labjournal/blob/master/labjournal/FastQC_results/Trashless_fqc_December/V300021405_L02_2_1_Adaptless_Trashless_fastqc.html),
[3](http://htmlpreview.github.io/?https://github.com/regnveig/labjournal/blob/master/labjournal/FastQC_results/Trashless_fqc_December/V300021405_L02_3_1_Adaptless_Trashless_fastqc.html),
[4](http://htmlpreview.github.io/?https://github.com/regnveig/labjournal/blob/master/labjournal/FastQC_results/Trashless_fqc_December/V300021405_L02_4_1_Adaptless_Trashless_fastqc.html).

## Выравнивание

```bash
boomer;
Logo "BWA Align";
THREADS=10;
NOSORT_FOLDER='/dev/datasets/FairWind/_results/December_Dinara/bam_nosort';
GENOME="/dev/datasets/FairWind/_db/hg19/hg19.fa";
mkdir -p $NOSORT_FOLDER;
r1=(/dev/datasets/FairWind/_results/December_Dinara/Trashless/*1_Adaptless_Trashless.fq.gz);
r2=(/dev/datasets/FairWind/_results/December_Dinara/Trashless/*2_Adaptless_Trashless.fq.gz);
for var in ${!r1[*]};
do {
start_time=$(StartTime);
fb=$(FileBase ${r1[var]});
out=$NOSORT_FOLDER/${fb::-22}_unsorted.bam;
bwa mem -t $THREADS -v 1 $GENOME ${r1[var]} ${r2[var]} | samtools view -O BAM -@ $THREADS - > $out;
echo "File "${fb::-22}" is done "$(Timestamp $start_time)"";
} done;
Seal $NOSORT_FOLDER
```

```bash
boomer;
Logo "Sort & Index"
THREADS=10;
OUTPUT_FOLDER='/dev/datasets/FairWind/_results/December_Dinara/bam_sorted';
mkdir -p $OUTPUT_FOLDER;
for file in /dev/datasets/FairWind/_results/December_Dinara/bam_nosort/*.bam;
do PicardCommandLine SortSam SO=coordinate I=$file O=$OUTPUT_FOLDER/$(FileBase $file).bam;
done;
BamIndex $OUTPUT_FOLDER/*.bam;
Seal $OUTPUT_FOLDER
```

# The 96

## Обрезка адаптеров

```bash
boomer;
Logo "Illuminaless (cutadapt)";
THREADS=10;
OUTPUT_FOLDER='/dev/datasets/FairWind/_results/December_Dinara/The_96/Illuminaless';
mkdir -p $OUTPUT_FOLDER;
illumina1='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC';
illumina2='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT';
r1=(/dev/datasets/ngs_data/Dinara_December/191201_M02435_0064_000000000-CGHCH/Data/Intensities/BaseCalls/*R1_001.fastq.gz);
r2=(/dev/datasets/ngs_data/Dinara_December/191201_M02435_0064_000000000-CGHCH/Data/Intensities/BaseCalls/*R2_001.fastq.gz);
for var in ${!r1[*]};
do {
start_time=$(StartTime);
r1_output=$OUTPUT_FOLDER/$(FileBase ${r1[var]})_Illuminaless.fq.gz;
r2_output=$OUTPUT_FOLDER/$(FileBase ${r2[var]})_Illuminaless.fq.gz;
fb=$(FileBase ${r1[var]});
report=$OUTPUT_FOLDER/${fb::-6}_report.txt;
cutadapt -m 8 -j $THREADS -a $illumina1 -A $illumina2 -o $r1_output -p $r2_output ${r1[var]} ${r2[var]} > $report;
echo "File "${fb::-6}" is done "$(Timestamp $start_time)"";
} done;
Seal $OUTPUT_FOLDER
```

## Выравнивание

```bash
boomer;
Logo "BWA Align";
THREADS=10;
NOSORT_FOLDER='/dev/datasets/FairWind/_results/December_Dinara/The_96/Bam_NoSort';
GENOME="/dev/datasets/FairWind/_db/hg19/hg19.fa";
mkdir -p $NOSORT_FOLDER;
r1=(/dev/datasets/FairWind/_results/December_Dinara/The_96/Illuminaless/*R1_001_Illuminaless.fq.gz);
r2=(/dev/datasets/FairWind/_results/December_Dinara/The_96/Illuminaless/*R2_001_Illuminaless.fq.gz);
for var in ${!r1[*]};
do {
start_time=$(StartTime);
fb=$(FileBase ${r1[var]});
out=$NOSORT_FOLDER/${fb::-20}_unsorted.bam;
bwa mem -t $THREADS -v 1 $GENOME ${r1[var]} ${r2[var]} | samtools view -O BAM -@ $THREADS - > $out;
echo "File "${fb::-20}" is done "$(Timestamp $start_time)"";
} done;
Seal $NOSORT_FOLDER
```

## Удаление дубликатов

```bash
boomer;
Logo "Duplicate Remove";
output_folder="/dev/datasets/FairWind/_results/December_Dinara/The_96/Bam_dupless";
mkdir -p $output_folder;
mkdir -p $output_folder/queryname_sorted;
for file in /dev/datasets/FairWind/_results/December_Dinara/The_96/Bam_NoSort/*.bam;
do PicardCommandLine SortSam SO=queryname I=$file O=$output_folder/queryname_sorted/$(FileBase $file).bam;
done;
Seal $output_folder/queryname_sorted;
mkdir -p $output_folder/queryname_dupless;
for file in $output_folder/queryname_sorted/*.bam;
do PicardCommandLine MarkDuplicates REMOVE_DUPLICATES=true M=$output_folder/queryname_dupless/$(FileBase $file)_metrics.txt I=$file O=$output_folder/queryname_dupless/$(FileBase $file).bam;
done;
Seal $output_folder/queryname_dupless;
mkdir -p $output_folder/coordinate_dupless;
for file in $output_folder/queryname_dupless/*.bam;
do PicardCommandLine SortSam SO=coordinate I=$file O=$output_folder/coordinate_dupless/$(FileBase $file)_sorted.bam;
done;
Seal $output_folder/coordinate_dupless;
```

## Variant Calling

```bash
boomer;
Logo "Variant Calling";
THREADS=10;
FA="/dev/datasets/FairWind/_db/hg19/hg19.fa";
OUTPUT_PATH="/dev/datasets/FairWind/_results/December_Dinara/The_96/VCF";
DEPTH=10;
QUALITY=30;
mkdir -p $OUTPUT_PATH;
for file in /dev/datasets/FairWind/_results/December_Dinara/The_96/Bam_dupless/coordinate_dupless/*.bam;
do bcftools mpileup --threads $THREADS -f $FA $file | bcftools call --threads $THREADS -cv -Ou | bcftools filter -i "DP>$DEPTH & %QUAL>$QUALITY" > $OUTPUT_PATH/$(FileBase $file)_DP"$DEPTH"_QUAL"$QUALITY".vcf;
done
```
