import gzip
import sys
import pandas as pd
import numpy as np
import string
import re
from multiprocessing import cpu_count, Pool
import functools

from blister import *

def Out(total_):
	print(f"Total: {total_}", end='\r')

def thread0(filename):

	seqs = pd.DataFrame(np.array([['-bridge-', 'GCTGAGG'], ['-egdirb-', 'CCTCAGC'], ['-gatc-', 'GATC']]), columns=['name', 'seq'])
	genome = '-genome-'

	number = 1000000
	lines = number * 4
	percent = number / 100
	
	print(f"Start: {filename}:", end="\n")
	
	input0 = blister_input_open(filename, 'rt')
	
	counter = 0
	total = 0
	tyk = 0
	main_list = []
	main_table = pd.DataFrame(columns=['count', 'mask'])

	for it in range(lines):
		
		line = input0.readline()
	
		counter += 1
		if counter == 5:
			counter = 1
	
		if (counter == 2):
			line = line[0:-1]
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
	
	shorted_list = list(set(main_list))
	for note in shorted_list:
		main_table = main_table.append(pd.Series([main_list.count(note), note], index=['count', 'mask']), ignore_index=True)
	
	main_table['count'] = main_table['count'].apply(lambda x: x / percent)
	
	table_bgg = main_table[main_table['mask'].str.contains("^[^o\[]*\-bridge\-\-gatc\-\-genome-.*$")]
	table_eg = main_table[main_table['mask'].str.contains("^[^o\[]*\-egdirb\-\-genome\-.*$")]
	
	print(f"Done: {filename}", end="\n")
	
	table_result = pd.DataFrame(columns=['filename', 'bgg', 'eg'])
	table_result = table_result.append(pd.Series([filename, table_bgg['count'].sum(), table_eg['count'].sum()], index=['filename', 'bgg', 'eg']), ignore_index=True)
	input0.close()
	return table_result

THREADS_NUM = cpu_count()
pool = Pool(THREADS_NUM)

filenames = blister_input(["/dev/datasets/FairWind/_results/cut/illuminaless/sample-*_R1_Illuminaless.fastq.gz"])
if not filenames: exit()
output_folder = blister_dir('/dev/datasets/FairWind/_results/iless_anal/')
if not output_folder: exit()

results = pool.map(thread0, filenames)

pool.close()
pool.join()
del pool

table_result = pd.DataFrame(columns=['filename', 'bgg', 'eg'])
table_result = table_result.append(results, ignore_index=True)

output_file = blister_output("./samples_1-1-9", output_folder, "statistics", "csv", rewrite=True)
table_result.to_csv(output_file, sep='\t', index=False)
