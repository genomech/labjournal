# Обработка экзомных данных

## Протокол

### Оценка качества входных данных (*FastQC*)

```bash
fastqc -o $output_dir -t $core_number $input_dir/*.fastq.gz
```
* Размер последовательностей: >=8 букв
* Сверхрепрезентированные последовательности: отсутствуют

### Обрезка адаптеров (*cutadapt*)

### Выравнивание и сортировка (*bwa*)

* Убрать неканонические хромосомы из bed-файла

```bash
bwa mem -t $core_number -v 1 $ref $input_r1 $input_r2 | samtools view -O BAM -@ $core_number - > $unsorted_bam;
```

### Удаление дубликатов (*PicardTools*, *Strandless*)

1. *PicardTools*: маркирование дубликатов.
2. *Strandless*: удаление разнонаправленных копий.

```bash
PicardCommandLine SortSam SO=queryname I=$unsorted_bam O=$temp_bam;
PicardCommandLine MarkDuplicates REMOVE_DUPLICATES=true M=$picard_metrics_txt I=$temp_bam O=$dupless_bam && rm -f $temp_bam;
python3 $strandless_dir/Strandless.py -f BAM -i $dupless_bam -o $strandless_bam -m $strandless_metrics_txt;
```

### Рекалибровка qual'ов (*GATK*)

Для обучения модели потребуются вариации в VCF формате (для человеческого генома - [dbSNP >132](https://ftp.ncbi.nih.gov/snp/organisms/)).
Контиги у dbSNP другие, поэтому потребуется лёгкий перепарсинг:

```python
# python3 script.py input.vcf.gz output.vcf.gz
import pandas as pd
import gzip
chunksize = 10 ** 7
input_filename, output_filename = argv[1], argv[2]
with gzip.open(input_filename, 'wt') as output_file, gzip.open(output_filename, 'rt') as input_file:
	x = True
	while x: 
		line = input_file.readline()
		output_file.write(line)
		x = line[0:2] == '##'
	for chunk in pd.read_csv(input_file, chunksize=chunksize, sep='\t', header=None):
		chunk[0] = chunk[0].apply(lambda x: "chr" + str(x))
		chunk = chunk[chunk[0] != 'chrMT']
		output_file.write(chunk.to_csv(index=False, sep='\t', header=None))
```

Затем файл нужно пережать bgzip'ом (это важно) и заиндексировать.
Не юзать *gatk IndexFeatureFile*, он забагованный и баг они до сих пор не пофиксили.
Монстр в строке - это грубая доработка, GATK не жрёт строки с N-ками в Ref и Alt.
Лучше вынести это в скрипт выше.

```bash
zcat $reparsed_vcf_gz | grep -P '^(#.*)|([^\t]*\t[^\t]*\t[^\t]*\t[^N\t]*\t[^N\t]\t.*)$' | bgzip -c > $bgzipped_vcf_gz;
tabix -p vcf $bgzipped_vcf_gz
```

Если в bam-файле нет `@RG', то их придётся создать:

```bash
PicardCommandLine AddOrReplaceReadGroups I=$strandless_bam O=$readgroups_bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
```

Рекалибровка (предварительно):

```bash
cd $gatk_dir;
./gatk BaseRecalibrator -I $readgroups_bam --known-sites $dbsnp_vcf -O $recalibrate_table -R $ref
```

### Просмотр bam-файлов глазом (*IGV*)

```bash
PicardCommandLine SortSam SO=coordinate I=$strandless_bam O=$final_bam;
samtools index $final_bam;
```

### Анализ bam-файлов (*BamQC*)

Подготовка capture для QualiMap:

```bash
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,0,"."}' $capture.bed > "$capture"_QualiMap.bed
```

* Доля выравненных последовательностей: >80%
* Доля дубликатов в экзонах: <30%
* Обогащение (среднее покрытие capture / среднее покрытие вне capture): >10
* Покрытие 95%: >=10
* Среднее покрытие (арифметическое): >= 70

Это сильно завышенные показатели, но к ним надо стремиться.

### 6. Идентификация вариантов (*GATK*, *freebayes*, *vcfallelicprimitives*).

* Left realignment: *freebayes* делает это по умолчанию, но авторы протокола советуют делать это перед *freebayes*. 
* По какой-то сраной причине *freebayes* бастует и не желает фильтровать варианты по QUAL'у.
Делать это приходится с помощью bcftools.
* *vcfallelicprimitives* разбивает многобуквенные замены.
Транс-варианты он пишет как 1/0 и 0/1.

**TODO:** Уточнить почему.

```bash
freebayes -0 --min-coverage $min_coverage --max-coverage $max_coverage -f $ref -t $exome_bed -b $final_bam | bcftools filter -i "QUAL > "$min_qual"" | vcflib vcfallelicprimitives > $vcf;
```

* Skip-limit (skip regions of extremely high coverage): стартовое значение 200 (уточнить по ходу работы).
* Минимальная доля альтернативных аллелей для гетерозигот: 0.25
* QUAL: начать с 18, подстраивать под контрольные метрики из предыдущих пунктов (конкретно: глубину покрытия в экзоме).

**TODO:** freebayes выдает значения QUAL, которые означают "надежность" найденного варианта. В файл пишутся все варианты, но их нужно фильтровать. Можно стартовать с QUAL >20, но надо уточнить, какие значения используют другие.

* DP для идеальных данных (см. контроль качества выше): >12-15.
* DP при среднем покрытим <70: >8-10.
* DP при анализе конкретного гена: не учитывать.

### 8. Аннотация вариантов (*Annovar*, *Ensembl VEP*)

```bash
perl $annovar_dir/table_annovar.pl $output_vcf $annovar_dir/humandb -buildver $genome_assembly -protocol knownGene,ensGene,refGene,abraom,AFR.sites.2015_08,ALL.sites.2015_08,AMR.sites.2015_08,ASN.sites.2012_04,avgwas_20150121,avsift,avsnp150,cadd13,cg69,clinvar_20190305,cosmic70,dann,dbnsfp35c,dbscsnv11,EAS.sites.2015_08,eigen,esp6500_all,EUR.sites.2015_08,exac03,fathmm,gene4denovo201907,gerp++,gme,gnomad211_genome,gwava,hrcr1,icgc21,intervar_20180118,kaviar_20150923,ljb26_all,mcap13,mitimpact24,MT_ensGene,nci60,popfreq_all_20150413,regsnpintron,revel,SAS.sites.2015_08,snp142 --operation g,g,g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f --remove --vcfinput --thread $threads;
```
Все описанные БД надо доставлять отдельно в анноваре.

### 9. Поиск по таблице

* Отбор по частоте в популяционных данных (берём максимальную из всех БД, но в процессе анализа смотрим данные по всем).
	* Если частое аутосомно-рецессивное заболевание: частота по БД <5%.
	* Если редкое: <0.05-0.5%. 
	* Нижняя граница частоты: 0.01%. При анализе всегда берётся чуть выше.

**TODO:** Нужно проверить, почему в наших таблицах не всегда совпадают частоты с указанными в БД, если это глюк *Ensembl VEP*, перейти на *Annovar*. Пример: PRODH (HGMD) - нет частоты. В других БД частота высокая, в 1000 геномов, например.

* Поиск вариантов
	* Поиск loss-of-function variants: frameshift, stop codon, splice-site.
	* Если ничего не нашлось, ищем missence.

### 10. Проверка найденного

* Проверка сырых данных в bam-файле (не последний ли/первый/альтернативный экзон).
* Перепроверка частот в популяции по разным БД.
* GnomAD/ExAC loss-of-function intolerant (pLI) score для гена. Сверять с мышиными фенотипами, литературой и проч.
* HGMD. Описана ли в этой позиции другая замена и какой её эффект.
* ClinVar. Ей доверия меньше, поэтому после остальных.
* OMIM. Сравнение известного фенотипа с фенотипом по БД.
* SIFT, PolyPhen и пр.: предсказание эффекта.

### 11. Прочие особенности анализа

* Если есть подозрение на ограниченный список генов, сначала анализируются варианты в этих генах.
* Если нет: полный анализ.

## Требуемая информация по мутации

* Референсный аллель
* Замена
* Плоидность
* QUAL (*freebayes*)
* Имя гена
* Импакт
* Максимальная частота встречаемости из всех БД, частоты по разным БД

**Примечание.** Отдельным пунктом надо отмечать присутствие мутации в разных БД и разброс по частотам.

* pLI-score
* OMIM
* HGMD
* ClinVar
* Номер экзона, общее число экзонов
* DP: суммарное + для каждого из вариантов, если гетерозигота
* Координата (сортировка)
* Все остальные данные, в т.ч. данные программ, предсказывающих эффекты сплайсинга, замен АК, etc.

## Ссылки

1. [Калькулятор патогенности вариантов нуклеотидной последовательности](http://calc.generesearch.ru/)
