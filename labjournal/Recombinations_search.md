# Поиск мутаций

Наши библиотеки были обработаны *bcftools* с параметрами Depth = 10, Qual = 30.

```
$ bcftools mpileup -f {genome} {input_filename} | bcftools call -cv -Ou | bcftools filter -i \"DP>{depth} & %QUAL>{quality}\" > {output_filename}
```

Файлы получились чересчур большими, поэтому было решено дополнительно отфильтровать их на экзом (несортированный bed-файл MedExome TargetCapture):

```
vcftools --vcf {input_filename} --bed {region} --out {output_filename} --recode
```

Далее библиотеки обрабатывались в VEP со всеми возможными параметрами.

## Нейрофиброматоз I типа у пациента 1-9

Гены (hg19):

| Gene | Chrom | Begin    | End      |
|:-----|:-----:|:--------:|:--------:|
| NF1  | chr17 | 29421945 | 29709134 |

### NF1

| Ref | Alt | Coord          |
|:---:|:---:|:---------------|
| T   | C   | chr17:29463271 |
| C   | T   | chr17:29472077 |
chr17:29553485 exon

exon 29623288
