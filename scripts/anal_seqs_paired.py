import sys
import pandas as pd
import numpy as np
import string
import re
from Bio import SeqIO

from blister import *

def thread0():

	num = 100000
	seqs = pd.DataFrame(np.array([['-bridge-', 'GCTGAGG'], ['-egdirb-', 'CCTCAGC'], ['-gatc-', 'GATC']]), columns=['name', 'seq'])
	genome = '-genome-'

	output_dir = blister_dir("/dev/datasets/FairWind/_results/paired", create=True)
	output_file0 = blister_output("paired_1-7.txt", output_dir, "stat", "csv", rewrite=True)

	main_list1 = []
	main_table = pd.DataFrame(columns=['count', 'mask'])
	total = 0

	with blister_read("/dev/datasets/ngs_data/ExoC_Belopuz/30-213832944/sample-1-7_R1_001.fastq.gz", mode='rt') as input_file:
		for record in SeqIO.parse(input_file, "fastq"):

			line = record.seq.__str__()

			for index, row in seqs.iterrows():
				line = line.replace(row['seq'], row['name'])
			for it in range(10):
				its = str(it + 1)
				line = re.sub("(^|-)[ATGCN]{" + its + "}($|-)", "--[" + its + "]--", line)
	
			line = re.sub(r"[ATGCN]{11,}", genome, line)
			# optimization
			line = line.replace('ome--gatc--gen', '')
			line = re.sub(r"^-gatc--genome-", genome, line)
			line = re.sub(r"-genome--gatc-$", genome, line)
	
			main_list1.append(line)
			total += 1
			if total == num: break

	main_list2 = []
	total = 0

	with blister_read("/dev/datasets/ngs_data/ExoC_Belopuz/30-213832944/sample-1-7_R2_001.fastq.gz", mode='rt') as input_file:
		for record in SeqIO.parse(input_file, "fastq"):

			line = record.seq.__str__()

			for index, row in seqs.iterrows():
				line = line.replace(row['seq'], row['name'])
			for it in range(10):
				its = str(it + 1)
				line = re.sub("(^|-)[ATGCN]{" + its + "}($|-)", "--[" + its + "]--", line)
	
			line = re.sub(r"[ATGCN]{11,}", genome, line)
			# optimization
			line = line.replace('ome--gatc--gen', '')
			line = re.sub(r"^-gatc--genome-", genome, line)
			line = re.sub(r"-genome--gatc-$", genome, line)
	
			main_list2.append(line)
			total += 1
			if total == num: break

	main_list_total = []

	with blister_timestamp("CONCAT") as start_time:
		for i in range(num):
			main_list_total += [main_list1[i] + " | " + main_list2[i]]
			blister_progressbar(i / num, start_time)
		blister_erase()
			
	with blister_timestamp("COUNT") as start_time:
		shorted_list = list(set(main_list_total))
		for note in shorted_list:
			main_table = main_table.append(pd.Series([main_list_total.count(note), note], index=['count', 'mask']), ignore_index=True)
		main_table['count'] = main_table['count'].apply(lambda x: x / (total / 100))
		main_table.sort_values(by=['count'], ascending=False).to_csv(output_file0, sep='\t', index=False)
		
thread0()
