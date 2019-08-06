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
Было требование добавить в файл записи из базы данных HGMD, генотип, bai (SAM index).
Она была выкачана с помощью [Ensembl BioMart](http://grch37.ensembl.org/biomart/martview/) (DATABASE: Ensembl Variation 97, DATASET: Human Short Variants, filter by Variant Source: HGMD-PUBLIC).

Вариант скрипта *vape.py* на 6.08.19:

```python
import pandas as pd
import re
from multiprocessing import cpu_count, Pool
import functools
import time
import pickle

def parse_header(info):
    re_list = re.split("\|", info)
    re_list[0] = re.findall('\S+$', re_list[0])[0]
    re_list[-1] = re.findall('^[^"]+', re_list[-1])[0]
    return re_list

def parse_line(line, vcf_header, csq_header):
    line_list = re.split("\t", line)
    info = line_list[7]
    sample = line_list[9]
    del line_list[slice(-3,-1)]
    del line_list[-1]
    line_list += re.split(":", sample)
    line_list[-1] = (line_list[-1])[0:-1]
    gt_list = re.split("/", line_list[-2])
    if len(gt_list) == 2:
        if gt_list[0] == gt_list[1]:
            line_list[-2] += " [homo]"
        else:
            line_list[-2] += " [hetero]"
    
    chrom, pos = line_list[0], line_list[1]

    csq_block = pd.DataFrame(columns=csq_header)
    csq_bigline = re.split("CSQ=", info)[1]
    csq_lines = re.split(",", csq_bigline)
    for csq_line in csq_lines:
        csq_line_list = re.split("\|", csq_line)
        if len(csq_line_list) == len(csq_header):
            csq_block = csq_block.append(pd.Series(csq_line_list, index=csq_header), ignore_index=True)
        else:
            csq_block = csq_block.append(pd.Series(['DAMAGED'] * len(csq_header), index=csq_header), ignore_index=True)
            print(' - DAMAGED string found:\nFull: %s\nPart: %s', end="\n")
        
    
    vcf_block = pd.DataFrame(columns=vcf_header)
    for it in range(len(csq_block.index)):
        vcf_block = vcf_block.append(pd.Series(line_list, index=vcf_header), ignore_index=True)
    
    table = pd.concat([vcf_block, csq_block], axis=1)
    return table

# __main__

THREADS_NUM = cpu_count()

in_path = '/dev/datasets/FairWind/dinara/vep/'
filename = '38_S4' # .vcf
out_path = './'
hgmd_path = './hgmd.db'

input_file = open(in_path + filename + ".vcf", 'r')
pool = Pool(THREADS_NUM)

print(f"Start file '{filename}' on {THREADS_NUM} threads...\nStarted: {time.ctime(time.time())}", end="\n")
start_time = time.time()

csq_header = list()
vcf_header = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'GT', 'PL']

while 1:
    line = input_file.readline()
    if (re.match("^#[^#].*$", line) != None):
        break
    if (re.match("^##INFO=<ID=CSQ.*$", line) != None):
        csq_header = parse_header(line)

main_table = pd.DataFrame(columns=(vcf_header + csq_header))

print(f"Headers are ready [%f sec]" % (time.time() - start_time), end="\n")
start_time = time.time()

results = pool.map(functools.partial(parse_line, vcf_header=vcf_header, csq_header=csq_header), input_file)
pool.close()
pool.join()

print(f"Parsing is done [%f sec]" % (time.time() - start_time), end="\n")
start_time = time.time()

total = len(results)
main_table = main_table.append(results, ignore_index=True)
del results

print(f"Merging {total} strings is done [%f sec]" % (time.time() - start_time), end="\n")
start_time = time.time()

with open(f"./{filename}.pickle", 'wb') as f:
    pickle.dump(main_table, f)

print(f"Pickle is done [%f sec]" % (time.time() - start_time), end="\n")
start_time = time.time()

hgmd_data = pd.read_csv(hgmd_path, sep=',')
hgmd_data.rename(columns={"Variant name":"HGMD", "Chromosome/scaffold name":"CHROM", "Chromosome/scaffold position start (bp)":"POS"}, inplace=True)
hgmd_data["CHROM"] = hgmd_data["CHROM"].apply(lambda x: "chr"+str(x))
hgmd_data.drop(columns={'Chromosome/scaffold position end (bp)', 'Variant source'}, axis=1, inplace=True)
hgmd_data["POS"] = hgmd_data["POS"].apply(lambda x: str(x))

main_table = pd.merge(hgmd_data, main_table, how='right', on=["CHROM","POS"])
base_new.sort_values(by=["CHROM","POS"], inplace=True)

print(f"Merging table with HGMD database [%f sec]" % (time.time() - start_time), end="\n")
start_time = time.time()

main_table.to_csv(out_path + filename + ".csv", sep='\t', index=False, mode='a')

print(f"Writing to file(s) '{out_path + filename}.csv' is done [%f sec]" % (time.time() - start_time), end="\n")
```
