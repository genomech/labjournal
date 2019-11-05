from lib.blister import *
import pandas as pd

input_files = Blister.Input(["/dev/datasets/FairWind/_results/check/*.txt"])
if not input_files: exit()

output_dir = Blister.Dir("/dev/datasets/FairWind/_results/check/")
if not output_dir: exit()

main_table = pd.DataFrame(index=list(range(13, 150)))

print(main_table)

for input_file in input_files:
	table = pd.read_csv(input_file, sep="\t")
	main_table[input_file] = table['count']

main_table.to_csv("/dev/datasets/FairWind/_results/check/full.txt")
