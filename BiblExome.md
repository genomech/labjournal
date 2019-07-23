# Обработка экзомной библиотеки

Лаборатории Воеводы потребовалось обработать большой объём входных данных.
Процедура во многом рутинная, но на всякий случай запишу все шаги.

1. Загрузка hg19, компиляция Bowtie 2 Index.

```
$ wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz
$ bowtie2-build --threads 10 ./hg19.fa ./hg19/hg19
```

2. Проверка данных с помощью *FastQC*.

```
$ fastqc -o [output_dir] -t 6 [input_dir]/*.fastq.gz
```

Результаты [здесь](./FastQC_results/fastqc_dinara_190722).

3. Выравнивание с помощью *bowtie2*, конвертирование в BAM и сортировка.

Предварительно:

```bash
#!/bin/bash

for var in 104_S3 111_S6 113_S2 117_S5 38_S4 98_S1 le1_S7 le2_S8 le3_S9 le4_S10 le5_S11 le6_S12
do
bowtie2 -x /dev/datasets/FairWind/hg19_small/hg19 -1 /dev/datasets/ngs_data/biblexome/"$var"_R1_001.fastq.gz \
-2 /dev/datasets/ngs_data/biblexome/"$var"_R2_001.fastq.gz --very-sensitive \
-S -p 10 | samtools view -bS -@ 10 -O BAM | samtools sort -@ 10 \
-O BAM - > /dev/datasets/FairWind/dinara/"$var".sorted.bam
echo Files "$var" are ready.
done
```
