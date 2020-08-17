# NIPT

## Статьи

1. [Noninvasive Prenatal Testing Using Next Generation Sequencing: Pilot Experience of the D.O. Ott Research Institute of Obstetrics, Gynecology and Reproductology](https://sci-hub.tw/https://link.springer.com/article/10.1134/S1022795419100053)
2. Киты для выделения cfDNA
	* MagMAX Cell-Free DNA Isolation Kit 
	* [Cell-Free DNA - QIAGEN](https://www.qiagen.com/us/products/discovery-and-translational-research/dna-rna-purification/dna-purification/cell-free-dna/)
3. [Noninvasive prenatal testing for fetal subchromosomal copy number variations and chromosomal aneuploidy by low‐pass whole‐genome sequencing](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6565572/)

## Протокол Ivaschenko et al.

### Приготовление

1. Кровь матери (венозная?)
2. Кровь поместить в пробирки с 0.5M EDTA, pH 8.0 (Greiner Bio One, Austria), перемешать.
Не позже чем через 4 ч после взятия образца центрифугировать на 2000g, 10 мин при 4°C.
Затем собрать плазму и центрифугировать при 16000g, 10 мин при 4°C.
3. Отобрать ДНК из плазмы с помощью MagMAX Cell-Free DNA Isolation Kit (ThermoFisher Scientific Inc., United States), следуя рекомендациям производителя.
Стадия включает end-repair, пришивание баркодов и амплификацию.
4. Приготовить библиотеки согласно протоколу Ion Plus Fragment Library Kit (Thermo FisherScientific Inc.) с некоторыми модификациями (какими?).
5. Измерить концентрацию ДНК с помощью Qubit dsDNA HSAssay Kit (Thermo Fisher Scientific Inc.) на флуориметре Qubit 2.0 (Invitrogen, United States).
6. Качество ДНК - 2200 TapeStation Instrument (Agilent Technologies, United States) для капиллярного электрофореза с помощью High Sensitivity D1K ScreenTape и High Sensitivity D1K Reagents (Agilent Technologies).
7. Приготовить библиотеки для секвенирования и загрузить на чипы, используя IonChef System (ThermoFisher Scientific Inc.) с Ion 540 Kit-Chef и Ion 540 Chips (ThermoFisher Scientific Inc.).
8. Секвенирование произведено на Ion Torrent S5 (ThermoFisher Scientific Inc.)

### Биоинформатика

1. Фильтрация ридов по длине и квалу (SAMtools).
2. Выравнивание
3. Определение пола плода
4. Фетальная фракция определена с помощью R-пакета SeqFF
5. Нормализация GC-контента
6. Определение анэуплоидии с помощью Z-score регрессии (R-пакет NIPTeR).
