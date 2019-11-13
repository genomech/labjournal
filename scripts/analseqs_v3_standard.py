import sys
import pandas as pd
import numpy as np
import string
import re
from Bio import SeqIO

from lib.blister import *

def thread0(filename, output_dir):

	seqs = pd.DataFrame(np.array([['-bridge-', 'GCTGAGG'], ['-egdirb-', 'CCTCAGC'], ['-gatc-', 'GATC']]), columns=['name', 'seq'])
	genome = '-genome-'

	output_file0 = Blister.Output(filename[1], output_dir, "stat", "csv", rewrite=True, index=filename[0])

	main_list = []
	main_table = pd.DataFrame(columns=['count', 'mask'])
	total = 0

	with Blister.Read(filename[1], 'rt', filename[0]) as input_file:
		with Blister.Timestamp(f"PARSING", filename[0]) as start_time:
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
	
	with Blister.Timestamp(f"MAKE TABLE", filename[0]) as start_time:
		shorted_list = list(set(main_list))
		for note in shorted_list:
			main_table = main_table.append(pd.Series([main_list.count(note), note], index=['count', 'mask']), ignore_index=True)
		main_table['count'] = main_table['count'].apply(lambda x: x / (total / 100))
		main_table.sort_values(by=['count'], ascending=False).to_csv(output_file0, sep='\t', index=False)

	with Blister.Timestamp(f"COUNT BGG/EG", filename[0]) as start_time:
		table_bgg = main_table[main_table['mask'].str.contains("^[^o\[]*\-bridge\-\-gatc\-\-genome-.*$")]
		table_eg = main_table[main_table['mask'].str.contains("^[^o\[]*\-egdirb\-\-genome\-.*$")]
		table_result = pd.DataFrame(columns=['filename', 'bgg', 'eg'])
		table_result = table_result.append(pd.Series([filename[1], table_bgg['count'].sum(), table_eg['count'].sum()], index=['filename', 'bgg', 'eg']), ignore_index=True)
	
	return table_result

Blister.Logo("Analysis of Sequences v3.0 [blister]")

input_filenames = Blister.Input(["/dev/datasets/ngs_data/ma_add/*.fastq.gz"])
if not input_filenames: exit()
output_dir = Blister.Dir("/dev/datasets/FairWind/_results/check", create=True)
if not output_dir: exit()
output_file = Blister.Output("/all", output_dir, "", "csv", rewrite=True)

with Blister.Threading(f"ANALYSIS") as pool:
	results = pool.map(functools.partial(thread0, output_dir=output_dir), enumerate(input_filenames))

with Blister.Timestamp(f"MERGE BGG/EG") as start_time:
	table_result = pd.DataFrame(columns=['filename', 'bgg', 'eg'])
	table_result = table_result.append(results, ignore_index=True)
	table_result.to_csv(output_file, sep='\t', index=False) 
