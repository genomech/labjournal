from lib.blister import *
import pandas as pd

input_files = Blister.Input(["/dev/datasets/FairWind/_results/November/adapters/*bridge.txt.csv"])
if not input_files: exit()

output_dir = Blister.Dir("/dev/datasets/FairWind/_results/November/adapters/")
if not output_dir: exit()

main_table = pd.DataFrame(index=list(range(8, 150)))

for input_file in input_files:
	table = pd.read_csv(input_file, sep="\t")
	main_table[input_file] = table['count']

main_table.to_csv("/dev/datasets/FairWind/_results/November/adapters/all_bridge.csv")
