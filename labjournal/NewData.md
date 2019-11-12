# Ноябрьские данные BGI HiC

## Образцы

Скачать можно [здесь](http://genedev.bionet.nsc.ru/site/hic_out/2019-11-09-BGI-ExoC/data/results/).

Команда:

```bash
wget -r ftp://genedev.bionet.nsc.ru/site/hic_out/2019-11-09-BGI-ExoC/data/results
```

| #   | Sample                    |
|:---:|:--------------------------|
| 2   | #8 (Нарышкина), DNase     |
| 4   | Паша                      |
| 5   | #3 (Васильева), DNase     |
| 6   | MNase 0.5 (длинные)       |
| 7   | #9 (Бикмухаметова), DNase | 
| 15  | MNase 2 (короткие)        |
| 19  | K562 – DNase              |

Адаптеры:

1. Bridge

```
 P        Biot
 |        |
 5-GCTGAGGGATC-3
3-TCGACTCC-5  
```

2. Blunt

```
5-CAGTGGCGAC-3
3-GTCACCG-5   
```

**TODO:**

* Количество ридов в каждом файле
* % адаптеров и расположение
* MNAse-Hi-C и DNAse-Hi-C прогнать через FastQC
* Обрезать адаптеры
* anal_seqs

## Первичный анализ

### Количество ридов в каждом файле

Команда:

```bash
for file in ./*/*.fq.gz; do ( length=$(echo "print("$(zcat $file | wc -l)" / 4)" | python3); echo "| "$file" | "$length" |"; ) done
```
| Библиотека                                 | Ридов     |
|--------------------------------------------|-----------|
| ./15/191107_X603_FCH5KNCCCX2_L5_15_1.fq.gz | 61403140  |
| ./15/191107_X603_FCH5KNCCCX2_L5_15_2.fq.gz | 61403140  |
| ./19/191107_X603_FCH5KNCCCX2_L5_19_1.fq.gz | 60389488  |
| ./19/191107_X603_FCH5KNCCCX2_L5_19_2.fq.gz | 60389488  |
| ./2/191107_X603_FCH5KNCCCX2_L5_2_1.fq.gz   | 59088991  |
| ./2/191107_X603_FCH5KNCCCX2_L5_2_2.fq.gz   | 59088991  |
| ./4/191107_X603_FCH5KNCCCX2_L5_4_1.fq.gz   | 25522262  |
| ./4/191107_X603_FCH5KNCCCX2_L5_4_2.fq.gz   | 25522262  |
| ./5/191107_X603_FCH5KNCCCX2_L5_5_1.fq.gz   | 138769927 |
| ./5/191107_X603_FCH5KNCCCX2_L5_5_2.fq.gz   | 138769927 |
| ./6/191107_X603_FCH5KNCCCX2_L5_6_1.fq.gz   | 15847105  |
| ./6/191107_X603_FCH5KNCCCX2_L5_6_2.fq.gz   | 15847105  |
| ./7/191107_X603_FCH5KNCCCX2_L5_7_1.fq.gz   | 121075216 |
| ./7/191107_X603_FCH5KNCCCX2_L5_7_2.fq.gz   | 121075216 |

### % адаптеров и расположение

Команда.
Нет, это не призыв сотоны из глубин ада, это всего лишь pipe в два параллельных процесса:

```bash
bridge="GCTGAGGGATC"; blunt="CAGTGGCGAC"; error=0.2; null="/dev/null"; output="/dev/datasets/FairWind/_results/November/adapters"; for file in ./*/*.fq.gz; do ( str=$(echo $file | tr '/' '_' ); str=${str:2:-6}; mkfifo pipe_$str; cat pipe_$str | (cutadapt -j 1 -b $bridge -O 9 -e $error -o $null - > $output/"$str"_bridge.txt) & (zcat $file | head -4000000) | tee pipe_$str | (cutadapt -j 1 -b $blunt -O 8 -e $error -o $null - > $output/"$str"_blunt.txt); rm pipe_$str; echo $file is done. ) & done
```

Подготовка таблиц:

```bash
for file in ./*bridge.txt; do (tail -146 $file > $file.csv) done
for file in ./*bridge.txt; do (tail -144 $file > $file.csv) done
```

![Bridge adapter](./scripts_results/November_cutadapt_reports_bridge.png)
![Blunt adapter](./scripts_results/November_cutadapt_reports_blunt.png)

Тот же график без выброса в 4 библиотеке:

![Blunt adapter without outlier](./scripts_results/November_cutadapt_reports_blunt_not4.png)

Таблица ODS, [если понадобится](./scripts_results/November_cutadapt_reports.ods).
