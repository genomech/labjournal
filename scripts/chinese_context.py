import gzip
import sys
import numpy as np
import string
from lib.blister import *
from Bio import SeqIO

Blister.Logo("Context Analyzator")

C_POSITIONS = 12
C_READ_LENGTH = 150
C_SEQ = "GCTGAGG"

filename = '/dev/datasets/FairWind/_results/November/Illuminaless/191107_X603_FCH5KNCCCX2_L5_2_1_Illuminaless.fq.gz'
length = len(C_SEQ)
poses = list(range(- C_POSITIONS, 0)) + list(range(length, length + C_POSITIONS))
nuc = list('ATGCN0')
a = np.zeros(shape=(C_READ_LENGTH, C_POSITIONS * 2, 6))

total = 0

with Blister.Read(filename, 'rt') as input_file, Blister.Timestamp(f"PARSING") as start_time:
	for record in SeqIO.parse(input_file, "fastq"):
		line = record.seq.__str__()
		tyk = line.find(C_SEQ)
		if tyk != -1:
			for i in range(len(poses)):
				pos_abs = tyk + poses[i]
				if (pos_abs < 0) or (pos_abs > 149) or (len(line) <= pos_abs): a[tyk, i, nuc.index('0')] += 1
				else: a[tyk, i, nuc.index(line[pos_abs])] += 1
		total += 1
		if total > 1000000: break

print(np.sum(a, axis = (1, 2)))
np.save('/dev/datasets/FairWind/_results/test.dump', a)
