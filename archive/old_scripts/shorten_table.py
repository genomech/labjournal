from lib.blister import *

import pandas as pd

def the_thread(block, output_dir):
	
	index = block[0]
	input_filename = block[1]
	
	output_filename = Blister.Output(input_filename, output_dir, "Shorten", "csv", rewrite=True)
	if not output_filename: return False
	
	with Blister.Timestamp("READ CSV", index) as start_time:
		table = pd.read_csv(input_filename, sep='\t', header=None, names=['chr', 'depth', 'reads', 'total', 'percent', 'trash0', 'trash1', 'trash2', 'trash3'])
	
	with Blister.Timestamp("PROCESSING & SAVE", index) as start_time:
		table = table[table['chr'] == "all"]
		table.drop(columns=['chr', 'reads', 'total', 'trash0', 'trash1', 'trash2', 'trash3'], axis=1, inplace=True)
		table.to_csv(output_filename, sep='\t', index=False)

Blister.Logo("Coverage Shorten")

input_filenames = Blister.Input(["/dev/datasets/FairWind/_results/38_S4--CoverageExome.txt"])
if not input_filenames: exit()

output_dir = Blister.Dir("/dev/datasets/FairWind/_results/", create=True)
if not output_dir: exit()

with Blister.Threading("COVERAGE") as pool:
	pool.map(functools.partial(the_thread, output_dir=output_dir), enumerate(input_filenames))
