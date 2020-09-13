# Анализ MG1-2

Последовательности:

| Pattern | Seq |
|:-----|:---------|
| a-bridge | aGCTGAGGgatc |
| egdirb-t | gatcCCTCAGCt |
| bridge | GCTGAGGgatc |
| egdirb | gatcCCTCAGC |
| a-blunt | aGCTGAGGgac |
| tnulb-t | gtcCCTCAGCt |
| blunt | GCTGAGGgac |
| egdirb | gtcCCTCAGC |
| illumina | AGATCGGAAG |

Команда:

```bash
./FastContext.py -p '{"a-bridge": "AGCTGAGGGATC", "egdirb-t": "GATCCCTCAGCT", "bridge": "GCTGAGGGATC", "egdirb": "GATCCCTCAGC", "a-blunt": "AGCTGAGGGAC", "tnulb-t": "GTCCCTCAGCT", "blunt": "GCTGAGGGAC", "egdirb": "GTCCCTCAGC", "illumina": "AGATCGGAAG"}' -i /dev/datasets/FairWind/_results/MG12_analysis/MG_2.fq.gz -o /dev/datasets/FairWind/_results/MG12_analysis/MG_2.csv 
```
