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

**Ген (hg19):** NF1, chr17:29,421,945-29,709,134

### Найденные перестройки

Результаты поиска [здесь](./scripts_results/1-9_NF1_search.csv).
Они были пропущены через VEP со всеми возможными параметрами.
Результаты [здесь](./scripts_results/1-9_NF1_vep.csv).
