# Обработка экзомной библиотеки

Лаборатории Воеводы потребовалось обработать большой объём входных данных.
Процедура во многом рутинная, но на всякий случай запишу все шаги.
За основу брались скрипты Полины.

1. Загрузка hg19, компиляция Bowtie 2 Index.

```
$ wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz
$ bowtie2-build --threads [core_number] ./hg19.fa ./hg19/hg19
```

2. Проверка данных с помощью *FastQC*.

```
$ fastqc -o [output_dir] -t [core_number] [input_dir]/*.fastq.gz
```

Результаты [здесь](./FastQC_results/fastqc_dinara_190722).

3. Выравнивание с помощью *bowtie2*, конвертирование в BAM и сортировка с помощью *samtools*.

```bash
#!/bin/bash

THREADS=[core_number]
B2_INDEX="[bowtie2_index]"
INPUT_PATH="[input_path]"
OUTPUT_PATH="[output_path_for_SAMs]"
SORTED_PATH="[output_path_for_sorted_BAMs]"
METRICS_PATH="[output_path_for_metrics]"
ERR_PATH="[output_path_for_error_logs]"

for var in 104_S3 111_S6 113_S2 117_S5 38_S4 98_S1 le1_S7 le2_S8 le3_S9 le4_S10 le5_S11 le6_S12
do
( bowtie2 -x $B2_INDEX -1 $INPUT_PATH/"$var"_R1_001.fastq.gz -2 $INPUT_PATH/"$var"_R2_001.fastq.gz \
        --very-sensitive -p $THREADS --met-file $METRICS_PATH/"$var"_metrics.txt | tee \
        $OUTPUT_PATH/"$var".sam | samtools view -bS -@ $THREADS - | samtools sort -@ $THREADS \
        -O BAM - > $SORTED_PATH/"$var".sorted.bam ) 3>&1 1>&2 2>&3 | tee $ERR_PATH/"$var"_stderr.log
echo Files "$var" are ready.
done
```

4. Быстрая проверка BAM-файлов.

```
$ samtools quickcheck -v -v -v ./*.bam > ./bad_bams.fofn
```

5. Компиляция индекса для bcftools из hg19.fa с помощью *faidx*:

```
$ samtools faidx ./hg19.fa
```

6. Наложение фильтров с помощью *bcftools*.

```bash
#!/bin/bash

THREADS=[core_number]
B2_INDEX="[bowtie2_index]"
OUTPUT_PATH="[output_path_for_VCFs]"
SORTED_PATH="[output_path_for_sorted_BAMs]"
DEPTH=20
QUALITY=30

for var in 104_S3 111_S6 113_S2 117_S5 38_S4 98_S1 le1_S7 le2_S8 le3_S9 le4_S10 le5_S11 le6_S12
do
bcftools mpileup --threads $THREADS -f $B2_INDEX $SORTED_PATH/"$var".sorted.bam | \
bcftools call --threads $THREADS -cv -Ou | bcftools filter -i "DP>$DEPTH & %QUAL>$QUALITY" \
> $OUTPUT_PATH/filterCalls_"$var"_DP"$DEPTH"_QUAL"$QUALITY"
echo File "$var" is ready.
done
```

7. Полученные файлы обрабатываются с помощью [VEP Ensembl, GRCH37](http://grch37.ensembl.org/Homo_sapiens/Tools/VEP), с подключением всех возможных инструментов.
Результат скачивается в формате *vcf*.

8. Дополнение полученных файлов.
Было требование добавить в файл записи из базы данных [HGMD](http://www.hgmd.cf.ac.uk/ac/index.php), генотип, bai (SAM index).

Предварительный вариант скрипта *vape.py*:

```python
import pandas as pd
import re

def parse_header(info):
    re_list = re.split("\|", info)
    re_list[0] = re.findall('\S+$', re_list[0])[0]
    re_list[-1] = re.findall('^[^"]+', re_list[-1])[0]
    return re_list

filename = './test_short.vcf'
input_file = open(filename, 'r')

csq_header = list()
vcf_header = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'FORMAT', 'SAMPLE']

while 1:
    line = input_file.readline()
    if (re.match("^#[^#].*$", line) != None): 
        break
    if (re.match("^##INFO=<ID=CSQ.*$", line) != None):
        csq_header = parse_header(line)

csq_table = pd.DataFrame(columns=csq_header)
vcf_table = pd.DataFrame(columns=vcf_header)

for line in input_file:
 
    line_list = re.split("\t", line)
    info = line_list[7]
    del line_list[7]

    # Разбираем INFO. Кроме CSQ, прочие значения Динаре вроде как не нужны, поэтому с ними не работаем
    
    csq_block = pd.DataFrame(columns=csq_header)
    csq_bigline = re.findall("[^=]+$", info)[0]
    csq_lines = re.split(",", csq_bigline)
    for csq_line in csq_lines:
        csq_line_list = re.split("\|", csq_line)
        csq_block = csq_block.append(pd.Series(csq_line_list, index=csq_header), ignore_index=True)
 
    csq_table = csq_table.append(csq_block, ignore_index=True)
    for it in range(len(csq_block.index)):
        vcf_table = vcf_table.append(pd.Series(line_list, index=vcf_header), ignore_index=True)

main_table = pd.concat([vcf_table, csq_table], axis=1)
main_table.to_excel("./excel.xlsx", index=False)
```
