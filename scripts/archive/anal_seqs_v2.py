import sys
import pandas as pd
import numpy as np
import string
import re
from Bio import SeqIO

from blister import *

def thread0(filename):

	seqs = pd.DataFrame(np.array([['-bridge-', 'GCTGAGG'], ['-egdirb-', 'CCTCAGC'], ['-gatc-', 'GATC']]), columns=['name', 'seq'])
	genome = '-genome-'

	output_dir = blister_dir("/dev/datasets/FairWind/_results/t78", create=True)
	output_file0 = blister_output(filename, output_dir, "stat", "csv", rewrite=True)

	main_list = []
	main_table = pd.DataFrame(columns=['count', 'mask'])
	total = 0

	with blister_read(filename, mode='rt') as input_file:
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
	
			main_list.append(line)
			total += 1
			if total == 1000000: break
	
	with blister_timestamp(f"Writing table {filename}") as start_time:
		shorted_list = list(set(main_list))
		for note in shorted_list:
			main_table = main_table.append(pd.Series([main_list.count(note), note], index=['count', 'mask']), ignore_index=True)
		main_table['count'] = main_table['count'].apply(lambda x: x / (total / 100))
		main_table.sort_values(by=['count'], ascending=False).to_csv(output_file0, sep='\t', index=False)

	with blister_timestamp(f"Counting bgg/eg {filename}") as start_time:
		table_bgg = main_table[main_table['mask'].str.contains("^[^o\[]*\-bridge\-\-gatc\-\-genome-.*$")]
		table_eg = main_table[main_table['mask'].str.contains("^[^o\[]*\-egdirb\-\-genome\-.*$")]
		table_result = pd.DataFrame(columns=['filename', 'bgg', 'eg'])
		table_result = table_result.append(pd.Series([filename, table_bgg['count'].sum(), table_eg['count'].sum()], index=['filename', 'bgg', 'eg']), ignore_index=True)
	
	return table_result

input_filenames = blister_input(["/dev/datasets/ngs_data/ExoC_Belopuz/30-213832944/sample-1-7_R*_001.fastq.gz", "/dev/datasets/ngs_data/ExoC_Belopuz/30-213832944/sample-1-8_R*_001.fastq.gz"])
if not input_filenames: exit()
output_dir = blister_dir("/dev/datasets/FairWind/_results/t78", create=True)
if not output_dir: exit()
output_file = blister_output("./all", output_dir, "", "csv", rewrite=True)

with blister_threading("Analysis of Sequences") as pool:
	results = pool.map(thread0, input_filenames)

table_result = pd.DataFrame(columns=['filename', 'bgg', 'eg'])
table_result = table_result.append(results, ignore_index=True)

table_result.to_csv(output_file, sep='\t', index=False)
