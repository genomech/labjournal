# Голый землекоп

## Задача

Мы пытаемся клонировать в плазмиду ген голого землекопа [PARP1](https://benchling.com/s/seq-euMSNt6Bv4O61ESmvlmA/edit), который связан с репарацией.
В ходе клонирования мы столкнулись с проблемой ПЦР-амплификации этого гена.
Некоторе фрагменты кодирующей последовательности амплифицируются нормально, а некоторые - совсем не получается наработать.

Геном голого землекопа собран не очень хорошо, поэтому возможно, что последовательность, по которой мы заказываем праймеры, просто неверная.
Тем более, что в разных БД мы нашли как минимум два разных варианта этого гена (правда, один встречается почти во всех БД, а второй - только в одной какой-то старой, скорее всего он неправильный).

Нужно качать несколько RNA-seq запусков, которых на землекопе много, и выровнять на фрагмент генома, предположительно содержащий ген PARP1, чтобы мы могли перепроверить, правильно ли этот ген аннотирован.

Актуальный геном голого землекопа: *HetGla_female_1.0/hetGla2*.
Последовательность ДНК можно скачать с USCS или NCBI.

PARP1 (предположительно): 
* Начало: JH602120:477,932-479,471
* Конец: JH602120:447,109-447,455

Cкачать всю эту область +/- 1 MB, потом выровнять на неё какие-нибудь транскриптомные данны разных типов клеток землекопа.
Ещё у меня просьба выровнять транскриптом на два варианта кДНК PARP1 из двух разных БД:

1. [ENSEMBL](https://benchling.com/s/seq-euMSNt6Bv4O61ESmvlmA)
2. [NCBI](https://benchling.com/s/seq-BCmTPIyyjbKSEgc6j0T3)

Использовать hisat (если ему не обязательно подавать аннотацию генов).
Если он без неё не может работать, то обычным bowtie (или bwa) можно, но только убрать ограничение на расстояние между выравниваниями двух ридов в паре (по дефолту у bowtie там что-то около 800 букв) - из-за сплайсинга это расстояние может быть очень большим.

## RNA-Seq

Найденные на NCBI:

| Specimen             | Tissue | SRA        |
|:---------------------|:-------|:-----------|
| BF1                  | brain  | SRS899007  |
| BF1                  | ovary  | SRS898991  |
| BF2                  | brain  | SRS899005  |
| BF2                  | ovary  | SRS898990  |
| BM1                  | testis | SRS898987  |
| BM1                  | brain  | SRS898999  |
| BM2                  | testis | SRS898980  |
| BM2                  | brain  | SRS898997  |
| E-MTAB-4550:Sample1  | liver  | ERS1090459 |
| E-MTAB-4550:Sample2  | liver  | ERS1090470 |
| E-MTAB-4550:Sample46 | liver  | ERS1090495 |
| E-MTAB-4550:Sample50 | liver  | ERS1090500 |
| SF-mix1              | ovary  | SRS898989  |
| SF-mix2              | ovary  | SRS898988  |
| SF1                  | brain  | SRS899003  |
| SF2                  | brain  | SRS899001  |
| SM1                  | testis | SRS898979  |
| SM1                  | brain  | SRS898995  |
| SM2                  | testis | SRS898978  |
| SM2                  | brain  | SRS898993  |

Использованы:

| Specimen             | Tissue | SRA        | SRR |
|:---------------------|:-------|:-----------|:----|
| BF1                  | brain  | SRS899007  | SRR1959124 |
| BM1                  | testis | SRS898987  | SRR1959204, SRR1959205 |
| E-MTAB-4550:Sample1  | liver  | ERS1090459 | ERR1331665 |

## Выравнивание

Предварительно:

```bash
for h2_index in Ensembl_CCDS NCBI_PARP2 JH602120_1MB;
do {
dir="/dev/datasets/FairWind/_results/SandPuppy/"$h2_index"";
mkdir -p "$dir";
index_path="/dev/datasets/FairWind/_db/hetGla2/"$h2_index"/"$h2_index"";
input_dir="/dev/datasets/ngs_data/HetGla_RNAseq";
for reads in ERR1331665 SRR1959124 SRR1959204 SRR1959205;
do {
hisat2 -x "$index_path" -1 ""$input_dir"/"$reads"_1.fastq.gz" -2 ""$input_dir"/"$reads"_2.fastq.gz" -p 10 --rg-id $reads --rg SM:sample"$reads" --rg LB:lib"$reads" --rg PL:ILLUMINA --rg PU:unit"$reads" | samtools view -O BAM > $dir/"$reads".bam;
} done;
} done
```
Мердж двух ранов одного образца:

```bash
samtools merge -O BAM SRR1959204-5_merged.bam SRR1959204.bam SRR1959205.bam
```

Sort:

```bash
for h2_index in Ensembl_CCDS NCBI_PARP2 JH602120_1MB;
do {
dir="/dev/datasets/FairWind/_results/SandPuppy/_SortedPuppy/"$h2_index"";
mkdir -p "$dir";
input_dir="/dev/datasets/FairWind/_results/SandPuppy/"$h2_index"";
for reads in ERR1331665 SRR1959124 SRR1959204-5_merged;
do {
PicardCommandLine SortSam SO=coordinate I="$input_dir"/"$reads".bam O="$dir"/"$reads".bam;
samtools index "$dir"/"$reads".bam;
} done;
} done
```

## Результат

[Здесь](./scripts_results/LightenedPuppy.zip).
