# NIPT

## Статьи

1. [Noninvasive Prenatal Testing Using Next Generation Sequencing: Pilot Experience of the D.O. Ott Research Institute of Obstetrics, Gynecology and Reproductology](https://sci-hub.tw/https://link.springer.com/article/10.1134/S1022795419100053)
2. Киты для выделения cfDNA
	* MagMAX Cell-Free DNA Isolation Kit 
	* [Cell-Free DNA - QIAGEN](https://www.qiagen.com/us/products/discovery-and-translational-research/dna-rna-purification/dna-purification/cell-free-dna/)
3. [Noninvasive prenatal testing for fetal subchromosomal copy number variations and chromosomal aneuploidy by low‐pass whole‐genome sequencing](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6565572/)

---

4. [A method for noninvasive detection of fetal large deletions/duplications by low coverage massively parallel sequencing](https://pubmed.ncbi.nlm.nih.gov/23592436/).
5. [Performance Evaluation of NIPT in Detection of Chromosomal Copy Number Variants Using Low-Coverage Whole-Genome Sequencing of Plasma DNA](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4945049/)
6. [WISECONDOR: detection of fetal aberrations from shallow sequencing maternal plasma based on a within-sample comparison scheme](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3950725/)
7. [Detection of fetal subchromosomal abnormalities by sequencing circulating cell-free DNA from maternal plasma ](https://pubmed.ncbi.nlm.nih.gov/25710461/)
8. [Limited Clinical Utility of Non-invasive Prenatal Testing for Subchromosomal Abnormalities](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4716686/)

---

9. [Non-invasive prenatal testing of fetal whole chromosome aneuploidy by massively parallel sequencing](https://pubmed.ncbi.nlm.nih.gov/23299662/)
10. [High resolution global chromosomal aberrations from spontaneous miscarriages revealed by low coverage whole genome sequencing ](https://pubmed.ncbi.nlm.nih.gov/29525519/)
11. [Noninvasive prenatal testing of fetal aneuploidies by massively parallel sequencing in a prospective Chinese population ](https://pubmed.ncbi.nlm.nih.gov/23703459/)


## Протокол [Ivaschenko et al](https://sci-hub.tw/https://link.springer.com/article/10.1134/S1022795419100053).

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

## Протокол [Yu et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6565572/) -- микроперестройки

### Приготовление

1. 10 мл периферической крови, ЭДТА-содержащие пробирки (Sekisui, Tokyo, Japan).
2. Плазму отделить не позже чем через 4 ч после взятия образца (метод описан [здесь](https://pubmed.ncbi.nlm.nih.gov/30131598/)).
3. Хранить при -80°C.
4. Геномная ДНК из амниоцитов разрезана на 250bp, а cfDNA приготовлена с помощью MagMAX™ Cell‐Free DNA Isolation Kit (Applied Biosystems™ cat.: A29319).
5. Сборку библиотек, контроль качества и секвенирование - как в 9, 10, 11 (2.5 нг cfDNA или фрагментированной DNA, затем к ним были прикреплены баркодированные адаптеры (8 букв) и ПЦР).
6. Секвенирование на NextSeq 550AR (Annoroad Gene Technology Co., Ltd China)

### Биоинформатика

1. Требования: для амниоцитов - 7.5 млн ридов, длина 40bp, в перерасчёте на 0.1х генома (?). Для плазмы - 4.2 миллионов ридов, 40bp длиной, Q30 >95%
2. Выравнивание - BWA, удаление дубликатов (`... unique reads ...`)
3. Геном разделён на кусочки по 100Kb, риды внутри каждого кусочка подсчитаны.
4. Число ридов скорректировано в каждом окне с помощью LOWESS модели, учитывая GC-контент (дальше следует матан).
5. Проверка Z-scores на группе из 1000 образцов с низкой вероятностью хромосомных перестроек.
6. Добавочное преобразование Z-score с учётом данных "условно здоровых" образцов.

## Протокол [Yin et al.](https://www.nature.com/articles/s10038-018-0489-9) -- de novo SNV

-- 30-я неделя (надо бы пораньше, концентрация cfDNA будет меньше)

### Приготовление

1. 5 мл материнской периферической крови отцентрифугированы 10 мин при 4°C на 1600g. Кровяные клетки были отделены и процентрифугированы на 2500g в течение 10 мин, плазма - 2500g, 10 минут.
Фракции клеток и плазмы немедленно помещены на -80°C до дальнейшей обработки.
2. ДНК из материнских и фетальных клеток извлечена с помощью Amp Genomic DNA Kit (TIANGEN, China).
1 мл плазмы экстрагирован с помощью QIAamp Blood Kit (Qiagen, Germany), 10 нг ДНК вытащены для приготовления библиотеки, измерены Qubit 3.0 (Thermo Fisher Scientific).
3. Клеточная ДНК фрагментирована Covaris S220 (Covaris, Woburn, MA, USA) согласно протоколу производителя.
Клеточная и плазменная ДНК использованы для приготовления библиотек (блант-эндинг, А-тейлинг, пришивание адаптеров).
4. 17 циклов ПЦР.
5. ДНК гибридизована с SeqCap EZ Probes oligo pool (Roche Nimblegen, USA), со специально подобранной таргетной панелью, для выявления известных вариантов, вызывающих заболевания.
6. Библиотеки оценены с помощью Agilent 2100 Bioanalyzer (Agilent, USA) and by quantitative PCR.
7. 3 библиотеки с разными баркодами были секвенированы (2 х 100bp) на Illumina HiSeq 2500 sequencer (Illumina, USA) по инструкциям производителя.

**Почему три?**

### Биоинформатика

1. Фильтрация ридов с более чем 5% N и QUAL30 <50%
2. BWA (`aln -o1-e63 -i15 -L -l31 -k1 -t6`), удаление дубликатов.
3. Псевдотетраплоидное генотипирование (“AAAA”, “AAAB”, “ABAA”, “ABAB”, “ABBB”, “BBAB”, “BBBB”) -- MMFF. Из них 4 (“AAAB”, “ABAA”, “ABBB”, and “BBAB”) используются для оценки доли фетальной ДНК (FC).
4. Матан.
5. Сэнгер.
