# Анализ MG1-2

Последовательности:

| Pattern | Seq |
|:-----|:---------|
| a-bridge-gatc | AGCTGAGGGATC |
| gatc-egdirb-t | GATCCCTCAGCT |
| bridge-gatc | GCTGAGGGATC |
| gatc-egdirb | GATCCCTCAGC |
| a-blunt | AGCTGAGGGAC |
| tnulb-t | GTCCCTCAGCT |
| blunt | GCTGAGGGAC |
| tnulb | GTCCCTCAGC |
| a-bridge | AGCTGAGG |
| egdirb-t | CCTCAGCT |
| bridge | GCTGAGG |
| egdirb | CCTCAGC |
| illumina | AGATCGGAAG |

Команда:

```bash
./FastContext.py -k 13 -p '{"a-bridge-gatc": "AGCTGAGGGATC", "gatc-egdirb-t": "GATCCCTCAGCT", "bridge-gatc": "GCTGAGGGATC", "gatc-egdirb": "GATCCCTCAGC", "a-blunt": "AGCTGAGGGAC", "tnulb-t": "GTCCCTCAGCT", "blunt": "GCTGAGGGAC", "tnulb": "GTCCCTCAGC", "a-bridge": "AGCTGAGG", "egdirb-t": "CCTCAGCT", "bridge": "GCTGAGG", "egdirb": "CCTCAGC", "illumina": "AGATCGGAAG"}' -i /dev/datasets/FairWind/_results/MG12_analysis/MG_2.fq.gz -o /dev/datasets/FairWind/_results/MG12_analysis/MG_2.csv 
```
