import pysam
import pandas as pd

samfile = pysam.AlignmentFile("/dev/datasets/FairWind/_results/cut/NEW/dupless/191107_X603_FCH5KNCCCX2_L5_5_dupless.bam", "r")

for read in samfile:
	read = read.tostring()
	lst = read.split('\t')
	unmapped = lst[5] == "*"
	if unmapped:
		print(read)

