# Обработка экзомной библиотеки

Лаборатории Воеводы потребовалось обработать большой объём входных данных.
Процедура во многом рутинная, но на всякий случай запишу все шаги.

1. Загрузка hg19, компиляция Bowtie 2 Index.

```
$ wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz
$ bowtie2-build --large-index --threads 10 ./hg19.fa ./hg19/hg19
```

2. Проверка данных с помощью *FastQC*.

```
$ fastqc -o [output_dir] -t 6 [input_dir]/*.fastq.gz
```

Результаты [здесь](./FastQC_results/fastqc_dinara_190722).

3. Выравнивание с помощью *bowtie2*.

Предварительно:

```
bowtie2 -x /dev/datasets/FairWind/hg19.fa -1 /dev/datasets/ngs_data/biblexome/104_S3_R1_001.fastq.gz \
-2 /dev/datasets/ngs_data/biblexome/104_S3_R2_001.fastq.gz --very-sensitive \
-S  /dev/datasets/FairWind/dinara/104_S3.gz -p 10 | samtools view -bS - | samtools sort - > 104_S3.sorted.bam
```
