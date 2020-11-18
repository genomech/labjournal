# Barracuda

## Сравнение Cutadapt и BWA

Cutadapt: 

| Metrics | Sequence Not Found | F Only Found | R Only Found | Both Found |
|:---|:---|:---|:---|:---|
| % | 7,75 | 18,08 | 13,11 | 61,04 |

BWA: 

| Metrics | Pair Unmapped | Pair Misplaced | F Mapped Properly, Mate Unmapped | R Mapped Properly, Mate Unmapped | F Mapped Properly, Mate Misplaced | R Mapped Properly, Mate Misplaced | Both Mapped Properly |
|:---|:---|:---|:---|:---|:---|:---|:---|
| % | 4,42 | 4,75 | 8,28 | 4,88 | 7,25 | 3,06 | 67,36 |

## Пары топовых баркодов

![Top_Mate_Rate_500](../data/graphs/Barracuda/Top_Mate_Rate_500.svg)

![Top-10_Mates_Rate_100](../data/graphs/Barracuda/Top-10_Mates_Rate_100.svg)

![Top-10_Mates_Absolute_100](../data/graphs/Barracuda/Top-10_Mates_Absolute_100.svg)

## Кластерный анализ первых 23 баркодов (R)

Красный - топовый баркод, имеющий самый больший процент сходства с остальными.
Жёлтый - те риды, которые отличаются от него на 1 букву, зелёный - на 2 буквы.
Различия проверялись алгоритмом Левенштейна.

Топовый баркод почему-то имеет все A на местах случайных букв о_О:

```
AActcAAgaaAAcgtAActaAAtcgAAgtaAAgt
```

![Cluster](../data/graphs/Barracuda/cluster_analysis.svg)

Было выяснено, что, несмотря на преобладание А в топовых баркодах, она наполовину меньше представлена на вырожденных позициях:

| |  |
|:---:|:---:|
| ![FreqPosF](../data/graphs/Barracuda/FreqsPosF.svg) | ![FreqPosF](../data/graphs/Barracuda/FreqsPosR.svg) |

Также оказалось, что распределение ридов по числу букв А отличается от ожидаемого:

| |  |
|:---:|:---:|
| ![NumberA_F](../data/graphs/Barracuda/NumberA_F.svg) | ![NumberA_R](../data/graphs/Barracuda/NumberA_R.svg) |
| ![NumberG_F](../data/graphs/Barracuda/NumberG_F.svg) | ![NumberG_R](../data/graphs/Barracuda/NumberG_R.svg) |
