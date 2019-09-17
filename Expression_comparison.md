# Сравнение экспрессии

Образцы коллег из Германии.

**TODO:** Прочитать про:

* hisat2 ([github repo](https://github.com/DaehwanKimLab/hisat2))
* cufflinks ([docs](https://cole-trapnell-lab.github.io/cufflinks/))
* CummeRbund

Установка:

1. hisat2 и cufflinks устанавливаются из стандартных репозиториев убунты.
2. cummeRbund устанавливатся с помощью R:

```R
install.packages("BiocManager")
BiocManager::install("cummeRbund")
```

Библиотеки были проверены с помощью *FastQC* ([данные здесь](http://htmlpreview.github.io/?https://github.com/regnveig/labjournal/blob/master/FastQC_results/FatimaQadri_analysis_190917.html)).
Было обнаружено незначительное количество К-меров (менее 0,1%), а также повторяющиеся буквы в начале ридов (первые 8-10 букв).
