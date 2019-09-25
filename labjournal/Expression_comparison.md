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

# Fatimannusa Qadri

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

## План

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

| Cuff ID     | NCBI ID      | WT  | Mutant     | Fold change | q        |
|:------------|:-------------|:---:|:----------:|:-----------:|:--------:|
| XLOC_014279 | NM_053621    | 0.0 | 317.915009 | inf         | 0.000914 |
| XLOC_014365 | -            | 0.0 | 3.170720   | inf         | 0.006213 |
| XLOC_019088 | -            | 0.0 | 3.067730   | inf         | 0.009537 |
| XLOC_010932 | -            | 0.0 | 2.055790   | inf         | 0.000914 |
| XLOC_011264 | NM_001127652 | 0.0 | 1.710050   | inf         | 0.000914 |
| XLOC_022004 | -            | 0.0 | 1.450680   | inf         | 0.000914 |
| XLOC_014084 | NM_031577    | 0.0 | 1.351500   | inf         | 0.000914 |
| XLOC_003081 | NM_019174    | 0.0 | 1.238400   | inf         | 0.000914 |
| XLOC_019549 | NM_001013247 | 0.0 | 1.127030   | inf         | 0.000914 |
| XLOC_013848 | NM_019183    | 0.0 | 1.071980   | inf         | 0.000914 |
| XLOC_011263 | NM_001127652 | 0.0 | 1.009560   | inf         | 0.029060 |
| XLOC_016051 | -            | 0.0 | 0.888906   | inf         | 0.000914 |
| XLOC_011265 | NM_053721    | 0.0 | 0.868083   | inf         | 0.000914 |
| XLOC_001184 | NM_001106361 | 0.0 | 0.830649   | inf         | 0.000914 |
| XLOC_007381 | -            | 0.0 | 0.815152   | inf         | 0.011067 |
| XLOC_005610 | -            | 0.0 | 0.756804   | inf         | 0.000914 |
| XLOC_004482 | -            | 0.0 | 0.683683   | inf         | 0.000914 |
| XLOC_012266 | -            | 0.0 | 0.597145   | inf         | 0.000914 |
| XLOC_006229 | NM_019344    | 0.0 | 0.491085   | inf         | 0.000914 |
| XLOC_005609 | -            | 0.0 | 0.478454   | inf         | 0.000914 |

| Cuff ID     | NCBI ID      | WT       | Mut | Fold change | q        |
|:------------|:-------------|:--------:|:---:|:-----------:|:--------:|
| XLOC_021385 | NM_022946    | 1.287030 | 0.0 | -inf        | 0.000914 |
| XLOC_008696 | -            | 1.202040 | 0.0 | -inf        | 0.000914 |
| XLOC_014363 | NM_031651    | 1.055120 | 0.0 | -inf        | 0.000914 |
| XLOC_010132 | -            | 1.019710 | 0.0 | -inf        | 0.000914 |
| XLOC_011045 | -            | 0.709557 | 0.0 | -inf        | 0.002437 |
| XLOC_009417 | NM_053789    | 0.674614 | 0.0 | -inf        | 0.000914 |
| XLOC_015564 | -            | 0.491164 | 0.0 | -inf        | 0.000914 |

Гены с изменением экспрессии в 2^(1..8) раз: [уменьшение](./scripts_results/genes_minus_diff.csv) и [увеличение](./scripts_results/genes_plus_diff.csv) экспрессии (CSV).

Гены, у которых исходная и/или конечная экспрессия больше 1000:

| Cuff ID     | NCBI ID      | WT          | Mutant     | Fold change | q       |
|:------------|:-------------|:-----------:|:-----------:|:---------:|:--------:|
| XLOC_013923 | NM_017016    | 9.969800    | 1020.400024 | 6.677350  | 0.000914 |
| XLOC_013547 | NM_013157    | 258.191986  | 1477.660034 | 2.516800  | 0.000914 |
| XLOC_006778 | NM_013144    | 334.798004  | 1469.750000 | 2.134210  | 0.000914 |
| XLOC_015152 | NM_021261    | 945.752991  | 3454.659912 | 1.869010  | 0.000914 |
| XLOC_020821 | NM_017072    | 424.191010  | 1410.849976 | 1.733780  | 0.006804 |
| XLOC_017303 | NM_182474    | 418.816986  | 1045.040039 | 1.319170  | 0.013562 |
| XLOC_005197 | NM_017233    | 493.239990  | 1131.560059 | 1.197950  | 0.006213 |
| XLOC_015733 | NM_053288    | 442.756989  | 1007.700012 | 1.186480  | 0.000914 |
| XLOC_004892 | -            | 1439.410034 | 2815.189941 | 0.967753  | 0.000914 |
| XLOC_003558 | -            | 1647.420044 | 2795.000000 | 0.762638  | 0.030886 |
| XLOC_008458 | NM_001079898 | 1767.699951 | 1081.689941 | -0.708582 | 0.047844 |
| XLOC_011710 | NM_053587    | 2532.290039 | 1504.050049 | -0.751593 | 0.049207 |
| XLOC_008451 | NM_173329    | 1837.780029 | 1089.339966 | -0.754512 | 0.031166 |
| XLOC_006340 | NM_017170    | 1115.989990 | 533.198975  | -1.065580 | 0.000914 |
| XLOC_005382 | NM_013105    | 1243.369995 | 578.583008  | -1.103660 | 0.002437 |
| XLOC_008791 | -            | 1077.089966 | 493.636993  | -1.125620 | 0.000914 |
| XLOC_014542 | NM_012556    | 2533.929932 | 878.020996  | -1.529050 | 0.000914 |

**TODO:** совместить гены крысёныша и CuffLinks ID. Сделано, результаты обновлены.

Для совмещения таблиц я впервые в жизни собственноручно пользовался поиском по бинарному древу:

```python
from lib.blister import *
import intervals as I

class Node: 
	def __init__(self, key): 
		self.left = None
		self.right = None
		self.val = key

def insert(root, node):
	if root is None: root = node 
	else:
		if root.val['interval'].lower < node.val['interval'].lower: 
			if root.right is None: root.right = node 
			else: insert(root.right, node) 
		else: 
			if root.left is None: root.left = node 
			else: insert(root.left, node)

def search(root, value): 
	if root is None: return None
	else:
		if root.val['interval'].is_connected(value): return root.val['ncbi_id']
		else:
			if root.val['interval'].lower < value.lower: return search(root.right, value) 
			else: return search(root.left, value)

def inorder(root): 
	if root: 
		inorder(root.left) 
		print(root.val)
		inorder(root.right) 

Blister.Logo("Sayeeda's Gene Finder")

with Blister.Timestamp("LOAD TRANSCRIPTOME") as start_time:
	transcriptome = pd.read_csv("/dev/datasets/FairWind/_results/Fatima/cuffmerge/merged.gtf", sep='\t', header=None, names=['chrom', 'hz1', 'hz2', 'begin', 'end', 'hz3', 'hz4', 'transcript_id', 'gene_id'])
	transcriptome.drop(columns=['hz1', 'hz2', 'hz3', 'hz4'], inplace=True, axis=0)
	transcriptome['transcript_id'] = transcriptome['gene_id'].apply(lambda x: x.split(" ")[3][1:-2])
	transcriptome['gene_id'] = transcriptome['gene_id'].apply(lambda x: x.split(" ")[1][1:-2])
	for col in ['begin', 'end']:
		transcriptome[col] = pd.to_numeric(transcriptome[col], downcast='integer', errors='raise')
	
	length = transcriptome.shape[0]
	intervals = []
	for it in transcriptome.iterrows():
		intervals += [I.IntInterval.closed(it[1]['begin'], it[1]['end'])]
		assert not I.IntInterval.closed(it[1]['begin'], it[1]['end']).empty
		Blister.ProgressBar(it[0] / length, start_time)
		
	transcriptome['interval'] = intervals
	transcriptome.drop(columns=['begin', 'end'], inplace=True, axis=0)
	del intervals
	Blister.Erase()

with Blister.Timestamp("LOAD GENOME") as start_time:
	genome = pd.read_csv("/dev/datasets/FairWind/_db/rn6.bed", sep='\t', header=None, names=['chrom', 'begin', 'end', 'ncbi_id', 'hz1', 'hz2', 'hz3', 'hz4', 'hz5', 'hz6', 'hz7', 'hz8'])
	genome.drop(columns=['hz1', 'hz2', 'hz3', 'hz4', 'hz5', 'hz6', 'hz7', 'hz8'], inplace=True, axis=0)
	for col in ['begin', 'end']:
		genome[col] = pd.to_numeric(genome[col], downcast='integer', errors='raise')

	chroms = list(set(genome['chrom'].to_list()))
	sorted_genome = dict()
	done = 0
	for chrom in chroms:
		table = genome[genome['chrom'] == chrom]
		trigger = False
		sorted_genome[chrom] = None
		for it in table.iterrows():
			it[1]['interval'] = I.IntInterval.closed(it[1]['begin'], it[1]['end'])
			it[1].drop(labels=['begin', 'end', 'chrom'], inplace=True)
			if trigger: insert(sorted_genome[chrom], Node(it[1]))
			else:
				sorted_genome[chrom] = Node(it[1])
				trigger = True
		del table
		done += 1
		Blister.ProgressBar(done / len(chroms), start_time)
	Blister.Erase()

with Blister.Timestamp("FIND GENES") as start_time:
	length = transcriptome.shape[0]
	ncbi = []
	for trans in transcriptome.iterrows():
		if trans[1]['chrom'] in sorted_genome.keys():
			ncbi += [search(sorted_genome[trans[1]['chrom']], trans[1]['interval'])]
			Blister.ProgressBar(trans[0] / length, start_time)
		else:
			ncbi += [None]
	transcriptome['ncbi_id'] = ncbi
	Blister.Erase()

with Blister.Timestamp("PROCESS TABLE & SAVE") as start_time:
	transcriptome = transcriptome['ncbi_id'].groupby(transcriptome['gene_id'])
	transcriptome = transcriptome.apply(set).apply(list).apply(lambda x: [i for i in x if i]).apply(lambda x: x[0] if x else None)
	transcriptome.to_csv("/dev/datasets/FairWind/_results/Fatima/gene_id.csv", sep='\t', index=True)
```

## Выводы

Общие впечатления -- что-то с иммунным ответом.

* В 2^6 раз увеличилась экспрессия гистидиндекарбоксилазы (NM_017016) -- единственный фермент, синтезирующий гистамин.
* С нуля до 317 стартанула мембраноассоциированная гуанилаткиназа. Есть данные, что она участвует в метаболизме пуринов и противовирусном ответе (?).
* Регулятор G-белка -- логично связывается с уровнем гуанилаткиназы (GDP-GTP).
* Упал до нуля DLG-ассоциированный белок, но возрос альфа-1 актин -- цитоскелет.
* В 2^7 раз возрос фактор роста фибробластов (NM_130754).
* В 2^6 раз возрос протоонкоген (NM_001256509)
* В 2^6 раз восзос тропонин (NM_017184). Какого хрена мускульный актин и тропонин вообще делают в печени?
* Упал до нуля ИЛ-17B -- ИЛ противовирусного ответа.
* Задействовано много мембранных белков (solute career)
* В 2^5 раз возрос фактор активации транскрипции (NM_012912)

# Andre Felipe Rodrigues

Аналогичная методика.

## Результаты

![График полной экспрессии](./scripts_results/Andre_Expression_all.svg)
![График полной экспрессии, log](./scripts_results/Andre_Expression_all_log.svg)
![График экспрессии q < 0.05](./scripts_results/Andre_Expression_all_q005.svg)
![График экспрессии q < 0.05, log](./scripts_results/Andre_Expression_all_q005_log.svg)

Гены, прекратившие работать:

| Gene ID     | NCBI ID      |       Day |   Night |  log2 (fold change) |   q       | Gene description |
|:------------|:-------------|:---------:|:-------:|:-------------------:|:---------:|:---|
| XLOC_005100 | -            | 189.105   |       0 |                -inf | 0.0232269 | - |
| XLOC_002266 | NM_001277154 | 114.4     |       0 |                -inf | 0.0232269 | RN serine/arginine repetitive matrix 2 (Srrm2) |
| XLOC_011074 | NM_024349    |  92.7481  |       0 |                -inf | 0.0232269 | RN adenylate kinase 1 (Ak1) |
| XLOC_017099 | NM_001013165 |  42.6763  |       0 |                -inf | 0.0232269 | RN leucine rich repeat containing 23 (Lrrc23) |
| XLOC_002641 | NM_080583    |  37.359   |       0 |                -inf | 0.0232269 | RN adaptor related protein complex 2 subunit beta 1 (Ap2b1) |
| XLOC_017100 | NM_001013165 |  37.0375  |       0 |                -inf | 0.0232269 | RN leucine rich repeat containing 23 (Lrrc23) |
| XLOC_005723 | -            |  16.8746  |       0 |                -inf | 0.0232269 | - |
| XLOC_005619 | NM_001271241 |  15.1175  |       0 |                -inf | 0.0232269 | RN heterogeneous nuclear ribonucleoproteins C1/C2-like (LOC100911576) |
| XLOC_017752 | NM_001025728 |  13.3092  |       0 |                -inf | 0.0232269 | RN SWI/SNF related, matrix associated, actin dependent regulator of chromatin, subfamily b, member 1 (Smarcb1) |
| XLOC_000883 | NM_001107463 |  11.3395  |       0 |                -inf | 0.0232269 | RN VPS37C, ESCRT-I subunit (Vps37c) |
| XLOC_004838 | -            |  11.0581  |       0 |                -inf | 0.0232269 | - |
| XLOC_002640 | -            |  11.0254  |       0 |                -inf | 0.0232269 | - |
| XLOC_005461 | -            |   5.85781 |       0 |                -inf | 0.0232269 | - |
| XLOC_006483 | NM_001106055 |   4.15659 |       0 |                -inf | 0.0232269 | RN MYC binding protein 2 (Mycbp2) |

Гены, начавшие работать:

| Gene ID     | NCBI ID      |   Day |    Night |  log2 (fold change) |   q       | Gene description |
|:------------|:-------------|:-----:|:--------:|:-------------------:|:---------:|:---|
| XLOC_000169 | NM_001304816 |     0 | 72.053   |                 inf | 0.0232269 | RN paternally expressed 3 (Peg3) |
| XLOC_000668 | NM_012700    |     0 | 50.3297  |                 inf | 0.0232269 | RN syntaxin 1B (Stx1b) |
| XLOC_013338 | NM_001109616 |     0 | 35.5273  |                 inf | 0.0232269 | RN family with sequence similarity 219, member A (Fam219a) |
| XLOC_017751 | NM_001025728 |     0 | 17.1773  |                 inf | 0.0232269 | RN SWI/SNF related, matrix associated, actin dependent regulator of chromatin, subfamily b, member 1 (Smarcb1) |
| XLOC_012607 | NM_001134716 |     0 | 17.0699  |                 inf | 0.0232269 | RN C-type lectin domain family 12, member A (Clec12a) |
| XLOC_010242 | -            |     0 | 14.8634  |                 inf | 0.0232269 | - |
| XLOC_015082 | NM_053331    |     0 | 13.4997  |                 inf | 0.0232269 | RN thioredoxin 2 (Txn2) |
| XLOC_017361 | NM_001107920 |     0 | 12.8896  |                 inf | 0.0232269 | RN mitogen activated protein kinase kinase kinase 7 (Map3k7) |
| XLOC_006079 | NM_001107262 |     0 | 12.5023  |                 inf | 0.0232269 | RN spalt-like transcription factor 2 (Sall2) |
| XLOC_016089 | NM_001105997 |     0 | 12.4545  |                 inf | 0.0232269 | RN mitochondrial ribosomal protein L1 (Mrpl1) |
| XLOC_002121 | -            |     0 |  4.16469 |                 inf | 0.0232269 | - |
| XLOC_015909 | NM_012794    |     0 |  2.79158 |                 inf | 0.0232269 | RN glycosylation dependent cell adhesion molecule 1 (Glycam1) |
| XLOC_016036 | NM_001108137 |     0 |  2.75054 |                 inf | 0.0232269 | RN membrane frizzled-related protein (Mfrp) |

**Примечание.**
Транскрипты NM_001025728 находятся в обеих таблицах (изменение сплайсинга?).

Гены со значительным (2^(1..8)) изменением экспрессии:

1. Увеличение:

| Gene ID     | NCBI ID      |   Day    |    Night |  log2 (fold change) |   q       | Gene description |
|:------------|:-------------|:--------:|:--------:|:-------------------:|:---------:|:---|
| XLOC_007781 | NM_001101000 | 28.6775  | 569.142  |             4.3108  | 0.0232269 | RN microtubule-associated protein, RP/EB family, member 2 (Mapre2) |
| XLOC_017757 | NM_031138    |  2.30196 |  19.4957 |             3.08222 | 0.0232269 | RN ubiquitin-conjugating enzyme E2B (Ube2b) |
| XLOC_015224 | NM_022215    | 19.7019  |  64.983  |             1.72173 | 0.0232269 | RN glycerol-3-phosphate dehydrogenase 1 (Gpd1) |
| XLOC_017755 | -            | 89.2291  | 215.293  |             1.27071 | 0.0232269 | - |
| XLOC_014627 | NM_053595    | 13.0296  |  28.4926 |             1.1288  | 0.0232269 | RN placental growth factor (Pgf) |
| XLOC_001103 | NM_001193569 | 35.525   |  76.4545 |             1.10577 | 0.0232269 | RN serum/glucocorticoid regulated kinase 1 (Sgk1) |
| XLOC_001784 | NM_031834    |  9.21026 |  19.2096 |             1.06052 | 0.0232269 | RN sulfotransferase family 1A member 1 (Sult1a1) |

2. Уменьшение:

| Gene ID     | NCBI ID   |   Day    |    Night  |  log2 (fold change) |   q       | Gene description |
|:------------|:----------|:--------:|:---------:|:-------------------:|:---------:|:---|
| XLOC_006726 | -         | 422.419  | 114.129   |            -1.88801 | 0.0232269 | - |
| XLOC_006086 | NM_031056 |  22.7404 |   9.93753 |            -1.1943  | 0.0232269 | RN matrix metallopeptidase 14 (Mmp14) |
| XLOC_009955 | NM_030832 | 348.916  | 153.316   |            -1.18637 | 0.0232269 | RN fatty acid binding protein 7 (Fabp7) |

## Сравнение по зверям

Полные данные [здесь](./scripts_results/andre_s_beasts.csv).

Выбраны следующие кандидаты:

1. Sgk1 (XLOC_001103, NM_001193569) -- наиболее вероятный кандидат:

![Sgk1 expression](./scripts_results/andre_Sgk1.svg)

2. Mmp14 (XLOC_006086, NM_031056):

![Mmp14 expression](./scripts_results/andre_Mmp14.svg)

3. Sult1a1 (XLOC_001784, NM_031834) -- наименее вероятный кандидат, но тем не менее тоже может быть задействован:

![Sult1a1 expression](./scripts_results/andre_Sult1a1.svg)

**Примечание.**
Восьмой зверь выбивается по двум первым генам.
Я думаю, это какой-то неправильный зверь.
