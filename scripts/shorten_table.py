from lib.blister import *

import pandas as pd

def the_thread(block, output_dir):
	
	index = block[0]
	input_filename = block[1]
	
	output_filename = blister_output(input_filename, output_dir, "Shorten", "csv", rewrite=True)
	if not output_filename: return False
	
	with blister_timestamp("READ CSV", index) as start_time:
		table = pd.read_csv(input_filename, sep='\t', header=None, names=['chr', 'depth', 'reads', 'total', 'percent', 'trash0', 'trash1', 'trash2', 'trash3'])
	
	with blister_timestamp("PROCESSING & SAVE", index) as start_time:
		table = table[table['chr'] == "all"]
		table.drop(columns=['chr', 'reads', 'total', 'trash0', 'trash1', 'trash2', 'trash3'], axis=1, inplace=True)
		table.to_csv(output_filename, sep='\t', index=False)

blister_logo("Coverage Shorten 60M")

input_filenames = blister_input(["/dev/datasets/FairWind/_results/20-120M/coverage/*.txt"])
if not input_filenames: exit()

output_dir = blister_dir("/dev/datasets/FairWind/_results/20-120M/shorten_tables", create=True)
if not output_dir: exit()

with blister_threading("COVERAGE") as pool:
	results = pool.map(functools.partial(the_thread, output_dir=output_dir), enumerate(input_filenames))

print(results)
