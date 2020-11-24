# Анализ 20 образца

Последовательности:

| Pattern | Seq |
|:---|:---|
| linker | CGCGATATCTTATCTGAC |
| reknil | GTCAGATAAGATATCGCG |
| bridge-gatc  | GCTGAGGGATC |
| gatc-egdirb  | GATCCCTCAGC |
| blunt  | GCTGAGGGAC |
| tnulb  | GTCCCTCAGC |
| newblunt | CAGTGGCGAC |
| tnulbwen | GTCGCCACTG |
| bridge  | GCTGAGG |
| egdirb  | CCTCAGC |
| illumina  | AGATCGGAAG |

```bash
./FastContext.py -k 10 -p '{"linker": "CGCGATATCTTATCTGAC", "reknil": "GTCAGATAAGATATCGCG", "bridge-gatc": "GCTGAGGGATC", "gatc-egdirb": "GATCCCTCAGC", "blunt": "GCTGAGGGAC", "tnulb": "GTCCCTCAGC", "newblunt": "CAGTGGCGAC", "tnulbwen": "GTCGCCACTG", "bridge": "GCTGAGG", "egdirb": "CCTCAGC", "illumina": "AGATCGGAAG"}' -i /dev/datasets/FairWind/_results/AllBridgesArticle/s20_R1.fastq.gz -o /dev/datasets/FairWind/_results/AllBridgesArticle/s20_R1.csv
```

Контроль - экзом:

```bash
./FastContext.py -m 1000000 -k 10 -p '{"linker": "CGCGATATCTTATCTGAC", "reknil": "GTCAGATAAGATATCGCG", "bridge-gatc": "GCTGAGGGATC", "gatc-egdirb": "GATCCCTCAGC", "blunt": "GCTGAGGGAC", "tnulb": "GTCCCTCAGC", "newblunt": "CAGTGGCGAC", "tnulbwen": "GTCGCCACTG", "bridge": "GCTGAGG", "egdirb": "CCTCAGC", "illumina": "AGATCGGAAG"}' -i /dev/datasets/ngs_data/JulyThirteen_BGI/A/200714_X578_FCHCJYWCCX2_L2_13A-1.fq.gz -o /dev/datasets/FairWind/.cloud/core/labjournal/data/tables/s20_analysis_20201124/200714_X578_FCHCJYWCCX2_L2_13A-1.csv
```
