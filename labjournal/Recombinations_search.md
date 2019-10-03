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

**Гены (hg19):**

* NF1, chr17:29,421,945-29,709,134
* NF2, chr22:29,999,545-30,094,589

### NF1

Результаты поиска [здесь](./scripts_results/1-9_NF1_search.csv).
Они были пропущены через VEP со всеми возможными параметрами.
Результаты [здесь](./scripts_results/1-9_NF1_vep.csv).

Ничего значимого не найдено.

### NF2

Результаты поиска [здесь](./scripts_results/1-9_NF2_search.csv).
Они были пропущены через VEP со всеми возможными параметрами.
Результаты [здесь](./scripts_results/1-9_NF2_vep.csv).

Ничего значимого не найдено.

## Ген CFTR у пациента 1-1

CFTR: chr7:117,120,148- 117,307,162

### Результаты

Миссенс-мутация в chr7:117,267,786.
Глубина покрытия участка - 16, мутация видна в 4 ридах (25%), после эмпирического анализа на ПЦР-дубликаты - в 1 из 5.
Данные частоты по популяциям отсутствуют во всех БД.
Скорее всего, ошибка секвенирования.

![Mutation IGV](./scripts_results/1-1_CFTR_missense1.png)

Все прочие мутации имеют низкий импакт и/или высокую частоту встречаемости в популяции.
[Полные результаты](./scripts_results/1-1_CFTR_vep.csv).

## Братья с умственной отсталостью

Были смержены таблицы, а затем вручную отобраны гены, связанные с умственной отсталостью и агрессией.
Все они были вручную проверены в IGV.

| Chrom | Position  | Ref | Alt | Genotype | Type | Symbol | Gene Name_x | Phenotypes_x | Total AF |
|:------|:----------|:---:|:---:|:---------|:-----|:-------|:------------|:-------------|:--------:|
| chr10 | 283578    | G   | T   | 0/1 | missense | ZMYND11 | Zinc finger MYND domain-containing protein 11 | Mental retardation, autosomal dominant | 0.00016 |
| chr10 | 131639111 | C   | T   | 0/1 | missense | EBF3 | Early B-cell factor 3, Microtubule-associated protein, RP/EB family, member 3 | Hypotonia, ataxia, and delayed development syndrome, Autosomal dominant | 0.00012 |
| chr16 | 3819294   | C   | T   | 0/1 | missense | CREBBP | CREB binding protein | Menke-Hennekam syndrome 1, 618332 (3); Rubinstein-Taybi syndrome 1, Autosomal dominant | 0.00823 |
| chr17 | 7750177   | TACCACCACCACCACCACCACCACCACCACCACCACCACC | TACCACCACCACCACCACCACCACCACCACCACCACCACCACCACC | 0/1 [hetero] | inframe insertion | KDM6B | Lysine-specific demethylase 6B | Neurodevelopmental disorder with coarse facies and mild distal skeletal abnormalities, Autosomal dominant | - |
| chr17 | 17697093  | CCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA | CCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA | 0/1 [hetero] | inframe deletion | RAI1 | Retinoic acid-induced gene 1 | Smith-Magenis syndrome, Autosomal dominant, Isolated cases | - |
| chr2 | 165984281  | T   | C   | 0/1 | missense | SCN3A | Sodium channel, voltage-gated, type III, alpha polypeptide | Epilepsy, familial focal, with variable foci,  Autosomal dominant; Epileptic encephalopathy, early infantile, Autosomal dominant | 0.00023 |
| chr20 | 57428804  | A   | G   | 0/1 | missense | GNAS | GNAS complex locus (guanine nucleotide binding protein (G protein), alpha stimulating activity polypeptide 1) | McCune-Albright syndrome, somatic, mosaic; Osseous heteroplasia, progressive, Autosomal dominant; Pituitary adenoma, multiple types, somatic; Pseudohypoparathyroidism, Autosomal dominant; | 0.00699 |
| chr22 | 51135989  | GTT | G,GCCCCTT,GCCCCGCGCCCGGCCCCTT | 1/1 | splice acceptor, frameshift | SHANK3 | SH3 and multiple ankyrin repeat domains 3 | Phelan-McDermid syndrome, Autosomal dominant; Schizophrenia, Autosomal dominant  | - |
