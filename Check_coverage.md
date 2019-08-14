# Проверка покрытия экзомных районов ридами библиотеки

## Основное

Выяснить, насколько каждая буква районов в bad track (Roch MedExome target capture) покрывается ридами наших библиотек.
Это нужно для корректировки методики секвенирования.

**TODO**: Поискать программы, которые вычисляют coverage.

## Ход работы

1. Target capture можно скачать [отсюда](https://sequencing.roche.com/en/products-solutions/by-category/target-enrichment/hybridization/seqcap-ez-medexomekit.html) (раздел Design files).
Были использованы hg19 capture targets.

2. Для анализа coverage была выбрана утилита *bedtools*.
Плюсы её в том, что она есть в репозиториях убунты.

```
$ sudo apt install bedtools
```

*bedtools* оказался весьма глючной тулзой.
В общем, описываю порядок колдунства:

* Вытащить из заголовка sam-файла порядок хромосом.
Что-то типа такого: `@SQ	SN:chr1	LN:249250621`.
Очистить их от мусора (оставив только имена) и сохранить в *names.txt*.

* Отсортировать bed-файл, используя *names.txt*:

```
$ bedtools sort -faidx ./names.txt -i ./MedExome_hg19_capture_targets.bed > MedExome_hg19_capture_targets.sorted.bed
```

* сделать из hg19 fa.fai-файл (подробно о компиляции индекса [здесь](./BiblExome.md)).
* подключить fa.fai и отсортированный bed-файл:

```
bedtools coverage -hist -sorted -g /dev/datasets/FairWind/hg19.fa.fai -a /dev/datasets/FairWind/trimmed/MedExome_hg19_capture_targets.sorted.bed -b /dev/datasets/FairWind/trimmed/bams/sample-1-1.sorted.bam > ./coverage_report.txt
```

Ты можешь попробовать посчитать coverage без -sorted, мой юный друг, но помни -- в этом случае bam-файл уходит в оперативку полностью и процесс рискует убицца.
Да, чуваки как-то не подумали.
