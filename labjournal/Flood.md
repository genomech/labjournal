# Поток

## Калибровка по данным Динары

**Образец:** 104_S3

**Read length:** 100bp

![img](./scripts_results/Dinara_Calibri_coverage_100.svg)

| Parameter       | 20000000 | 30000000 | 40000000 | 50000000 | 60000000 |
|:----------------|:---------|:---------|:---------|:---------|:---------|
| Non-coverage, % | 0.24     | 0.15     | 0.12     | 0.1      | 0.08     |
| Coverage 50%    | 32       | 47       | 63       | 79       | 95       |
| Coverage 75%    | 21       | 31       | 42       | 53       | 64       |
| Coverage 90%    | 13       | 20       | 28       | 35       | 43       |
| Coverage 95%    | 9        | 15       | 21       | 26       | 32       |
| Twilight coverage | 39.01  | 58.37    | 77.79    | 97.33    | 116.7    |


Approximately 90M *paired* reads for coverage 50, and 140M reads for coverage 80.

## New pipeline

### CNNScoreVariants

Стандартные настройки из коробки, образец `DCSAN1-QUAR1`:

| Type | Variants, % |
|:-----|:------------|
| PASS | 97.56 |
| SNP, tranche 99.90-99.95 | 0.41 |
| SNP, tranche 99.95-100.00 | 1.38 |
| INDEL, tranche 99.00-99.40 | 0.2 |
| INDEL, tranche 99.40-100.00 | 0.41 |
