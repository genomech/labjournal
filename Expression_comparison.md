# Сравнение экспрессии

Образцы коллег из Германии.

**TODO:** Прочитать про:

* hisat2 ([github repo](https://github.com/DaehwanKimLab/hisat2))
* cufflinks ([docs](https://cole-trapnell-lab.github.io/cufflinks/))
* CummeRbund

Установка:

1. hisat2 и cufflinks устанавливаются из стандартных репозиториев убунты.
2. cummeRbund устанавливатся с помощью R:

```R
install.packages("BiocManager")
BiocManager::install("cummeRbund")
```

Библиотеки были проверены с помощью *FastQC* ([данные здесь](http://htmlpreview.github.io/?https://github.com/regnveig/labjournal/blob/master/FastQC_results/FatimaQadri_analysis_190917.html)).
Было обнаружено незначительное количество К-меров (менее 0,1%), а также повторяющиеся буквы в начале ридов (первые 8-10 букв).

[Индекс hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) генома крысы. 
С этим индексом *hisat2* крашнулся с segfault, поэтому я взял с него подписку о невыезде до выяснения обстоятельств дела.
Решил скачать нативный геном *R. norvegicus* с UCSC ([ссылка](https://hgdownload.soe.ucsc.edu/goldenPath/rn6/bigZips/)) и скомпилить его ручками.

```bash
hisat2-build ./rn6.fa ./rn6/rn6
```

Новый индекс вроде работает.
Старый был обвинён в измене и очищен удалением во славу индейских богов.

Статья с протоколом: [Nature Protocols](https://www.nature.com/articles/nprot.2012.016).

## Предварительный план

1. Выравнивание

```bash
#!/bin/bash

THREADS=12
HISAT2_INDEX="/dev/datasets/FairWind/_db/rn6/rn6"
INPUT_PATH="/dev/datasets/ngs_data/SeqDataTransfer/FatimannusaQadri/"
SAM_PATH="/dev/datasets/FairWind/_results/Fatima/sam"
BAM_PATH="/dev/datasets/FairWind/_results/Fatima/bam"

mkdir -p $SAM_PATH
mkdir -p $BAM_PATH

for var in '001_S13' '002_S14' '003_S15' '004_S16' '005_S17' '006_S18' '007_S19' '008_S20' '009_S21' '010_S22'
do
hisat2 -x $HISAT2_INDEX \
	-1 $INPUT_PATH/MB_FQ_"$var"_R1_001.fastq.gz \
	-2 $INPUT_PATH/MB_FQ_"$var"_R2_001.fastq.gz \
	-p $THREADS | tee $SAM_PATH/MB_FQ_"$var".sam | samtools view -bS \
	-@ $THREADS -O BAM - > $BAM_PATH/MB_FQ_"$var".bam;
echo $var is ready.
done

md5sum $SAM_PATH/*.sam > $SAM_PATH/all.md5;
chmod 555 -R $SAM_PATH;
md5sum $BAM_PATH/*.bam > $BAM_PATH/all.md5;
chmod 555 -R $BAM_PATH;

echo Dirs are sealed.
```

Средний процент выравнивания - 97%.

2. *Cufflinks* привереда и требует сортированные библиотеки.
Досортировываем и запускаем cufflinks.

```bash
#!/bin/bash

THREADS=12
BAM_PATH="/dev/datasets/FairWind/_results/Fatima/bam"
SORTED_PATH="/dev/datasets/FairWind/_results/Fatima/sorted"
CUFFLINKS_PATH="/dev/datasets/FairWind/_results/Fatima/cufflinks"

# Досортировка

mkdir -p $SORTED_PATH

for var in '001_S13' '002_S14' '003_S15' '004_S16' '005_S17' '006_S18' '007_S19' '008_S20' '009_S21' '010_S22'
do
samtools sort -@ $THREADS -O BAM $BAM_PATH/MB_FQ_"$var".bam > $SORTED_PATH/MB_FQ_"$var"_sorted.bam;
echo $var is sorted.
done

md5sum $SORTED_PATH/*.bam > $SORTED_PATH/all.md5;
chmod 555 -R $SORTED_PATH;

echo Sorted are sealed.

# Cufflinks

mkdir -p $CUFFLINKS_PATH

for var in '001_S13' '002_S14' '003_S15' '004_S16' '005_S17' '006_S18' '007_S19' '008_S20' '009_S21' '010_S22'
do
mkdir -p $CUFFLINKS_PATH/$var;
cufflinks -p $THREADS --library-type fr-firststrand \
	-o $CUFFLINKS_PATH/$var \
	$SORTED_PATH/MB_FQ_"$var"_sorted.bam;
echo $var is cuffled.
done

md5deep -lr $CUFFLINKS_PATH/* > $CUFFLINKS_PATH/all.md5;
chmod 555 -R $CUFFLINKS_PATH;

echo Cufflinks are sealed.
```

3. Мёржим транскрипты.

```bash
ls -1 /dev/datasets/FairWind/_results/Fatima/cufflinks/**/transcripts.gtf > /dev/datasets/FairWind/_results/Fatima/cuffmerge/transcripts_list.txt
cuffmerge -p 12 \
	-s /dev/datasets/FairWind/_db/rn6.fa \
	-o /dev/datasets/FairWind/_results/Fatima/cuffmerge/ \
	/dev/datasets/FairWind/_results/Fatima/cuffmerge/transcripts_list.txt
```

4. Запускаем *cuffdiff*.

	*  **Library Type** подбирался специально для наших библиотек.
	*  **Dispersion Method** и **Library Normalization Method** идут дефолтные, но я на всякий случай вынес их в команду, если потребуется настройка.

```bash
mkdir -p /dev/datasets/FairWind/_results/Fatima/cuffdiff
cd /dev/datasets/FairWind/_results/Fatima/sorted
cuffdiff -p 12 --library-type fr-firststrand \
	--dispersion-method pooled --library-norm-method geometric \
	-o /dev/datasets/FairWind/_results/Fatima/cuffdiff \
	/dev/datasets/FairWind/_results/Fatima/cuffmerge/merged.gtf \
	MB_FQ_006_S18_sorted.bam,MB_FQ_007_S19_sorted.bam,MB_FQ_009_S21_sorted.bam,MB_FQ_010_S22_sorted.bam \
	MB_FQ_001_S13_sorted.bam,MB_FQ_002_S14_sorted.bam,MB_FQ_003_S15_sorted.bam,MB_FQ_004_S16_sorted.bam,MB_FQ_005_S17_sorted.bam,MB_FQ_008_S20_sorted.bam
```
