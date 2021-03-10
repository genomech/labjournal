from Bio import SeqIO
from Levenshtein import distance as d
import gzip
import pandas as pd
import numpy as np

Bridge = "CTCAGC"
Ranger = 0

with gzip.open("/dev/datasets/ngs_data/ExoC_Belopuz/30-213832944/sample-1-1_R1_001.fastq.gz", "rt") as handle:
	for record in SeqIO.parse(handle, "fastq"):
		lst = []
		for length in range(len(Bridge) - Ranger, len(Bridge) + 1 + Ranger):
			strings = [record.seq[item : item + length].__str__() for item in range(len(record.seq) - length + 1)]
			strings = [d(Bridge, item) for item in strings] + ([0] * (len(record.seq) - len(strings)))
			lst += [strings]
		lst = np.mean(np.transpose(np.array(lst)), axis=1)
		print(lst)
		break
