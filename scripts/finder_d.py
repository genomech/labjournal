from lib.blister import *
import pandas as pd

def the_thread(block, output_dir):
	index = block[0]
	input_filename = block[1]
	
	output_filename = Blister.Output(input_filename, output_dir, "Filtered", "csv", rewrite=True, index=index)
	if not output_filename: return
	
	with Blister.Timestamp("READ CSV", filename_1=input_filename, index=index) as start_time:
		data = pd.read_csv(input_filename, sep='\t', low_memory=False)
		data['gnomAD_AF'] = pd.to_numeric(data['gnomAD_AF'], errors='coerce')
	
	with Blister.Timestamp("FILTER TABLE", index=index) as start_time:
		data = data[((data['IMPACT'] == "HIGH") | (data['IMPACT'] == "MODERATE")) & ((data['gnomAD_AF'] < 0.1) | (data['gnomAD_AF'].isna()))]
		
		all_genes = set(data['Gene'].to_list())
		saved = []
		for gene in all_genes:
			#if len(set(data[data['Gene'] == gene]['POS'].to_list())) > 1:
			if len(set(data[data['Gene'] == gene]['Location'].to_list())) > 1:
				saved += [gene]
			
		saved = set(saved)
		
		data = data[data['Gene'].isin(saved)]
		
		#data.sort_values(by=["CHROM","POS"], inplace=True)
		data.sort_values(by=['Location'], inplace=True)
	
	with Blister.Timestamp("SAVE TABLE", filename_1=output_filename, index=index) as start_time:
		data.to_csv(output_filename, sep='\t', index=False)
		del data
		del all_genes
		del saved

Blister.Logo("Finder D")

input_filenames = Blister.Input(["/dev/datasets/FairWind/_results/dinara/exome_polina/*.txt"])
#input_filenames = Blister.Input(["/dev/datasets/FairWind/_results/dinara/nepokocannoe_me/*.csv"])
if not input_filenames: exit()

#output_dir = Blister.Dir("/dev/datasets/FairWind/_results/dinara/find_stuff_my/", create=True)
output_dir = Blister.Dir("/dev/datasets/FairWind/_results/dinara/find_stuff_polina/", create=True)
if not output_dir: exit()

with Blister.Threading("PROCESSING FILES", THREADS_NUM=4) as pool:
	pool.map(functools.partial(the_thread, output_dir=output_dir), enumerate(input_filenames))

Blister.Seal(output_dir)
