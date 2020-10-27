# FastContext

FastQ Context Analysis Tool.

## Requirements

1. Python 3.8+
2. Python packages:

	* BioPython
	* contextlib
	* multiprocessing
	* argparse
	* bz2, gzip
	* functools
	* json
	* pandas

## Usage

| Parameter | Comment |
|:----------|:--------|
| -i INPUT_FILENAME | Fastq input file (may be gzipped or bzipped) |
| -o OUTPUT_FILENAME | Output CSV file |
| -p PATTERNS | Patterns to look for, plain JSON format: `{"pattern_1": "seq_1", "pattern_2": "seq_2"}` |
| -k KMER_SIZE | Max unrecognized K-mer size. Default = 0 |
| -u UNRECOGNIZED | Long unrecognized sequences replacement. Default = `genome` |
| -m MAX_READS | Max reads number to analyze (0 - no limit). Default = 0 |
| -f RATE_FLOOR | Min rate to write pattern in the table. Default = 0.01 |
| -@ THREADS | Threads number. Default = `cpu_count()` |

Example:

```bash
./FastContext.py \
	-k 13 \
	-p '{"bridge": "GCTGAGG", "egdirb": "CCTCAGC", "blunt": "GCTGAGGGAC", "tnulb": "GTCCCTCAGC", "illumina": "AGATCGGAAG"}' \
	-i input_file.fastq.gz \
	-o output.csv 
```

