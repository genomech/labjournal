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
