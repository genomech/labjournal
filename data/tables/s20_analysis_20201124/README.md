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
./FastContext.py -k 10 -p '{"linker": "CGCGATATCTTATCTGAC", "reknil": "GTCAGATAAGATATCGCG", "bridge-gatc": "GCTGAGGGATC", "gatc-egdirb": "GATCCCTCAGC", "blunt": "GCTGAGGGAC", "tnulb": "GTCCCTCAGC", "newblunt": "CAGTGGCGAC", "tnulbwen": "GTCGCCACTG", "bridge": "GCTGAGG", "egdirb": "CCTCAGC", "illumina": "AGATCGGAAG"}' -i /dev/datasets/FairWind/_results/AllBridgesArticle/s20_R1.fastq.gz /dev/datasets/FairWind/_results/AllBridgesArticle/s20_R1.csv
```
