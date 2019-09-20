# Сравнение экспрессии

Образцы коллег из Германии.

**TODO:** Прочитать про:

* hisat2 ([github repo](https://github.com/DaehwanKimLab/hisat2))
* cufflinks ([docs](https://cole-trapnell-lab.github.io/cufflinks/))
* CummeRbund

Установка:

1. hisat2 и cufflinks устанавливаются из стандартных репозиториев убунты.
2. cummeRbund устанавливатся с помощью R (под sudo!):

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

## cummeRbund

Первое.
Не запечатывать папки до первого запуска *cummeRbund*, т.к. он создаёт SQLite DB в папке.

После долгих мучений был получен следующий код.
Он создаёт CSV из db.

```R
library(cummeRbund)
cuff <- readCufflinks("/dev/datasets/FairWind/_results/Fatima/cuffdiff")
gene_diff_data <- diffData(genes(cuff))
write.table(as.data.frame(gene_diff_data),file="/dev/datasets/FairWind/genes3.csv", quote=F,sep="\t",row.names=F)
```

## Промежуточные результаты

![График полной экспрессии](./scripts_results/Expression_all.svg)
![График полной экспрессии, log](./scripts_results/Expression_all_log.svg)
![График экспрессии q < 0.05](./scripts_results/Expression_all_q005.svg)
![График экспрессии q < 0.05, log](./scripts_results/Expression_all_q005_log.svg)

Было решено отобрать следующие гены:

* Топ по fold change (и логарифму абсолютной экспрессии);
* q < 0.05 с уровнем экспрессии выше 1000.

Гены, у которых экспрессия появилась или пропала:

| Gene ID     | WT       | Mutant   | Fold change | q        |
|:------------|:--------:|:--------:|:-----------:|:--------:|
| XLOC_014279 | 0.000000 | 317.9150 | inf         | 0.000914 |
| XLOC_014365 | 0.000000 | 3.170720 | inf         | 0.006213 |
| XLOC_019088 | 0.000000 | 3.067730 | inf         | 0.009537 |
| XLOC_010932 | 0.000000 | 2.055790 | inf         | 0.000914 |
| XLOC_011264 | 0.000000 | 1.710050 | inf         | 0.000914 |
| XLOC_022004 | 0.000000 | 1.450680 | inf         | 0.000914 |
| XLOC_014084 | 0.000000 | 1.351500 | inf         | 0.000914 |
| XLOC_003081 | 0.000000 | 1.238400 | inf         | 0.000914 |
| XLOC_019549 | 0.000000 | 1.127030 | inf         | 0.000914 |
| XLOC_013848 | 0.000000 | 1.071980 | inf         | 0.000914 |
| XLOC_011263 | 0.000000 | 1.009560 | inf         | 0.029060 |
| XLOC_016051 | 0.000000 | 0.888906 | inf         | 0.000914 |
| XLOC_011265 | 0.000000 | 0.868083 | inf         | 0.000914 |
| XLOC_001184 | 0.000000 | 0.830649 | inf         | 0.000914 |
| XLOC_007381 | 0.000000 | 0.815152 | inf         | 0.011067 |
| XLOC_005610 | 0.000000 | 0.756804 | inf         | 0.000914 |
| XLOC_004482 | 0.000000 | 0.683683 | inf         | 0.000914 |
| XLOC_012266 | 0.000000 | 0.597145 | inf         | 0.000914 |
| XLOC_006229 | 0.000000 | 0.491085 | inf         | 0.000914 |
| XLOC_005609 | 0.000000 | 0.478454 | inf         | 0.000914 |

| Gene ID     | WT       | Mutant   | Fold change | q        |
|:------------|:--------:|:--------:|:-----------:|:--------:|
| XLOC_008696 | 1.202040 | 0.000000 | -inf        | 0.000914 |
| XLOC_009417 | 0.674614 | 0.000000 | -inf        | 0.000914 |
| XLOC_010132 | 1.019710 | 0.000000 | -inf        | 0.000914 |
| XLOC_011045 | 0.709557 | 0.000000 | -inf        | 0.002437 |
| XLOC_014363 | 1.055120 | 0.000000 | -inf        | 0.000914 |
| XLOC_015564 | 0.491164 | 0.000000 | -inf        | 0.000914 |
| XLOC_021385 | 1.287030 | 0.000000 | -inf        | 0.000914 |

Гены с изменением экспрессии в 2^(1..8) раз: [уменьшение](./scripts_results/genes_minus_diff.csv) и [увеличение](./scripts_results/genes_plus_diff.csv) экспрессии (CSV).

Гены, у которых исходная и/или конечная экспрессия больше 1000:

| Gene ID     | WT          | Mutant     | Fold change | q value |
|:------------|:-----------:|:-----------:|:---------:|:--------:|
| XLOC_013923 | 9.969800    | 1020.400024 | 6.677350  | 0.000914 |
| XLOC_013547 | 258.191986  | 1477.660034 | 2.516800  | 0.000914 |
| XLOC_006778 | 334.798004  | 1469.750000 | 2.134210  | 0.000914 |
| XLOC_015152 | 945.752991  | 3454.659912 | 1.869010  | 0.000914 |
| XLOC_020821 | 424.191010  | 1410.849976 | 1.733780  | 0.006804 |
| XLOC_017303 | 418.816986  | 1045.040039 | 1.319170  | 0.013562 |
| XLOC_005197 | 493.239990  | 1131.560059 | 1.197950  | 0.006213 |
| XLOC_015733 | 442.756989  | 1007.700012 | 1.186480  | 0.000914 |
| XLOC_004892 | 1439.410034 | 2815.189941 | 0.967753  | 0.000914 |
| XLOC_003558 | 1647.420044 | 2795.000000 | 0.762638  | 0.030886 |
| XLOC_008458 | 1767.699951 | 1081.689941 | -0.708582 | 0.047844 |
| XLOC_011710 | 2532.290039 | 1504.050049 | -0.751593 | 0.049207 |
| XLOC_008451 | 1837.780029 | 1089.339966 | -0.754512 | 0.031166 |
| XLOC_006340 | 1115.989990 | 533.198975  | -1.065580 | 0.000914 |
| XLOC_005382 | 1243.369995 | 578.583008  | -1.103660 | 0.002437 |
| XLOC_008791 | 1077.089966 | 493.636993  | -1.125620 | 0.000914 |
| XLOC_014542 | 2533.929932 | 878.020996  | -1.529050 | 0.000914 |
