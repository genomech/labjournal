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

Вариант скрипта *vape.py* на 7.08.19 (надеюсь, что окончательный):

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
filenames = ['104_S3', '111_S6', '113_S2', '117_S5', '38_S4', '98_S1', 'le1_S7', 'le2_S8', 'le3_S9', 'le4_S10', 'le5_S11', 'le6_S12'] 
out_path = './'
hgmd_path = './hgmd.db'

print(f"\n=== VAPE 0.13 ===\n", end="\n")
start_time = time.time()

hgmd_data = pd.read_csv(hgmd_path, sep=',')
hgmd_data.rename(columns={"Variant name":"HGMD", "Chromosome/scaffold name":"CHROM", "Chromosome/scaffold position start (bp)":"POS"}, inplace=True)
hgmd_data["CHROM"] = hgmd_data["CHROM"].apply(lambda x: "chr"+str(x))
hgmd_data.drop(columns={'Chromosome/scaffold position end (bp)', 'Variant source'}, axis=1, inplace=True) 
hgmd_data["POS"] = hgmd_data["POS"].apply(lambda x: str(x))

print(f"HGMD table is done [%f sec]" % (time.time() - start_time), end="\n")

for filename in filenames:

    input_file = open(in_path + filename + ".vcf", 'r')
    pool = Pool(THREADS_NUM)
    
    print(f"\nStart file '{filename}' on {THREADS_NUM} threads...\nStarted: {time.ctime(time.time())}", end="\n")
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
    del pool
    
    print(f"Parsing is done [%f sec]" % (time.time() - start_time), end="\n")
    start_time = time.time()
    
    total = len(results)
    main_table = main_table.append(results, ignore_index=True)
    del results
    
    print(f"Merging {total} strings is done [%f sec]" % (time.time() - start_time), end="\n")
    start_time = time.time()

    main_table = pd.merge(hgmd_data, main_table, how='right', on=["CHROM","POS"])
    main_table.sort_values(by=["CHROM","POS"], inplace=True)
    
    print(f"Merging table with HGMD database [%f sec]" % (time.time() - start_time), end="\n")
    start_time = time.time()
    
    with open(f"./pickles/{filename}.pickle", 'wb') as f:
        pickle.dump(main_table, f)
    
    print(f"Pickling is done [%f sec]" % (time.time() - start_time), end="\n")
    start_time = time.time()
    
    main_table.to_csv(out_path + "csv/" + filename + ".csv", sep='\t', index=False, mode='w')
    
    print(f"Writing to CSV file is done [%f sec]" % (time.time() - start_time), end="\n")
    
    del csq_header
    del vcf_header
    del main_table
    
    input_file.close()

```

9. Индекс с BAM-файлов снимался простым скриптом:

```bash
#!/bin/bash

# bambai.sh
  
for INFILE in "$@"
do
echo Start $INFILE ...
samtools index $INFILE
echo Done.
done
```

```
$ ./bambai.sh ./*.bam
```

10. Последняя задача, которую поставила Динара -- поискать соответствия в нескольких библиотеках.
Не зря я сохранял пиклзы.

```python
import pandas as pd
import functools
import time
import pickle

in_path = './pickles/'
filenames = ['le5_S11', 'le6_S12'] # .pickle
out_path = './csv_merged/'

print(f"\n=== MARGO 0.1 ===\n", end="\n")
start_time = time.time()

with open(f"{in_path}{filenames[0]}.pickle", 'rb') as f1:
    table1 = pickle.load(f1)
    
print(f"Pickle {filenames[0]} is loaded [%f sec]" % (time.time() - start_time), end="\n")
start_time = time.time()

with open(f"{in_path}{filenames[1]}.pickle", 'rb') as f2:
    table2 = pickle.load(f2)
    
print(f"Pickle {filenames[1]} is loaded [%f sec]" % (time.time() - start_time), end="\n")
start_time = time.time()

table_merged = pd.merge(table1, table2, how='inner')

del table1
del table2

print(f"Merging tables is done [%f sec]" % (time.time() - start_time), end="\n")
start_time = time.time()

table_merged.to_csv(out_path + filenames[0] + "-" + filenames[1] + "_merged.csv", sep='\t', index=False, mode='w')
    
print(f"Writing to CSV file is done [%f sec]" % (time.time() - start_time), end="\n")

```

---

**P.S.** Если бы кто-то погуглил прежде чем кодить, то узнал бы, что есть великолепная [библиотека для парсинга VCF](https://pyvcf.readthedocs.io/en/latest/).
Достаточно шустрая и удобная.
Так что, дорогой читатель, если ты вдруг захочешь сделать то же, что и я -- используй её, а не предыдущие костыли.

## Уменьшение размера таблицы

В ходе работы возникла следующая проблема -- таблицы чересчур большие, их сложно открывать на слабых компах.
Было решено прочесать их на предмет повторяющихся транскриптов и оставить только самые длинные.

БД транскриптов была скачана с Биомарта.
Параметры: Ensembl Genes 97, Human genes (GRCh37.p13).
Атрибуты:

* Gene stable ID
* Gene stable ID version
* Transcript stable ID
* Transcript stable ID version
* Transcript start (bp)
* Transcript end (bp)

Подготовка БД производилась скриптом *trans_filter.py*:

```python
import pickle
import pandas as pd

with open("../_pickles/csv/Ensembl_transcripts_[all].pd.pickle", 'rb') as f:
    transcripts = pickle.load(f)
    
transcripts['Transcript end (bp)'] = transcripts['Transcript end (bp)'].apply(pd.to_numeric, errors='ignore')
transcripts['Transcript start (bp)'] = transcripts['Transcript start (bp)'].apply(pd.to_numeric, errors='ignore')

transcripts['length'] = transcripts['Transcript end (bp)'] - transcripts['Transcript start (bp)']

new_transcripts = pd.DataFrame(columns=transcripts.columns)

dick = transcripts['Gene stable ID'].unique()

did = 0
total = len(dick)

for cock in dick:
    table = transcripts[transcripts['Gene stable ID'] == cock]
    new_transcripts = new_transcripts.append(table.loc[table['length'].idxmax()], ignore_index=True)
    did += 1
    print("%.2f%%" % (did * 100 / total), end='\r')

print("", end='\n')

new_transcripts.to_csv("./Ensembl_maxlength.csv", sep="\t", index=False)

print("Done", end='\n')

```

Номера NCBI можно скачать [отсюда](http://genome.ucsc.edu/cgi-bin/hgTables).
Но не нужно.

Изменения вносились с помощью скрипта *trans_removal.py*:

```python
import pandas as pd
import time
import pickle

filenames = ['104_S3', '111_S6', '113_S2', '117_S5', '38_S4', '98_S1', 'le1_S7', 'le2_S8', 'le3_S9', 'le4_S10', 'le5_S11', 'le6_S12'] 
out_path = './'

# tables

start_time = time.time()

ensembl_all = pd.read_csv("/dev/datasets/FairWind/_db/Ensembl_maxlength.csv", sep='\t')
ensembl_max = ensembl_all

ensembl_all = ensembl_all.drop(['Gene stable ID version', 'Transcript stable ID', 'Transcript stable ID version', 'Transcript start (bp)', 'Transcript end (bp)', 'length'], axis=1)
ensembl_all.rename(columns={'Gene stable ID':'Gene'}, inplace=True)
ensembl_all['blacksign_inbase'] = True

ensembl_max = ensembl_max.drop(['Gene stable ID version', 'Transcript stable ID', 'Transcript start (bp)', 'Transcript end (bp)', 'length'], axis=1)
ensembl_max.rename(columns={'Gene stable ID':'Gene', 'Transcript stable ID version':'Feature'}, inplace=True)
ensembl_max['blacksign_max'] = True

print(f"Ensembl tables are done [%f sec]" % (time.time() - start_time), end="\n")

for filename in filenames:
	
	print(f"\nStart file {filename} ...", end='\n')
	
	# pickle
	
	start_time = time.time()
	
	with open(f"/dev/datasets/FairWind/_results/dinara/pickles/{filename}.pickle", 'rb') as f:
		main_table = pickle.load(f)
	
	order = main_table.columns.to_list()
	
	print(f"Big data is done [%f sec]" % (time.time() - start_time), end="\n")
	
	# merge
	
	start_time = time.time()
	
	main_table = pd.merge(ensembl_all, main_table, how='right', on=["Gene"])
	main_table = pd.merge(ensembl_max, main_table, how='right', on=["Gene", "Feature"])
	main_table[['blacksign_max', 'blacksign_inbase']] = main_table[['blacksign_max', 'blacksign_inbase']].fillna(False)
	main_table = main_table.loc[~((main_table['blacksign_max'] == False) & (main_table['blacksign_inbase'] == True))]
	main_table = main_table.drop(['blacksign_max', 'blacksign_inbase'], axis=1)
	main_table = main_table[order]
	main_table.sort_values(by=["CHROM","POS"], inplace=True)
	
	print(f"Merging is done [%f sec]" % (time.time() - start_time), end="\n")
	
	# pickling
	
	start_time = time.time()
	    
	with open(f"/dev/datasets/FairWind/_results/dinara/new_pickles/{filename}_light.pickle", 'wb') as f:
		pickle.dump(main_table, f)
	    
	print(f"Pickling is done [%f sec]" % (time.time() - start_time), end="\n")
	
	# write to file
	
	start_time = time.time()
	
	main_table.to_csv(f"/dev/datasets/FairWind/_results/dinara/csv_light/{filename}_light.csv", sep='\t', index=False, mode='w')
	del main_table
	
	print(f"Writing to file is done [%f sec]" % (time.time() - start_time), end="\n")
```

Слияние родственников производилось вышеописанным *margo.py*, без изменений.

## Фильтрация Борна

Потребовалось отфильтровать некоторые записи в таблицах.
Заболевание рецессивное, вероятнее всего в нашем случае это компаунд гетерозиготы. 
Исходим из условий, что скорее всего 2 различные замены в одном гене приводят к появлению патологии. 

Схема работы:

1. Отфильтровываем все варианты  с частотой менее 1% и новые  по gnomad
2. Среди них отфильтровать с импактом High и Moderate  (на самом деле интересуют миссенсы, нонсенсы и сайт сплайсинговые) 
3. отфильтровать с 2-мя разными вариантами и более в одном гене. 

Экзомы:

* Полина -- 3, 49, 84, 95, 128, 141
* мои -- 38, 98, 104, 111, 113, 117.
