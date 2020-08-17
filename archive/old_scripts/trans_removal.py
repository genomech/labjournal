import pandas as pd
import time
import pickle

filenames = ['104_S3', '111_S6', '113_S2', '117_S5', '38_S4', '98_S1', 'le1_S7', 'le2_S8', 'le3_S9', 'le4_S10', 'le5_S11', 'le6_S12'] 
out_path = './'

# tables

start_time = time.time()

ensembl_all = pd.read_csv("/dev/datasets/FairWind/_db/Ensembl_maxlength.csv", sep='\t')
ensembl_max = ensembl_all

ensembl_all = ensembl_all.drop(['Gene stable ID version', 'Transcript stable ID', 'Transcript stable ID version', 'Transcript start (bp)', 'Transcript end (bp)', 'length'], axis=1)
ensembl_all.rename(columns={'Gene stable ID':'Gene'}, inplace=True)
ensembl_all['blacksign_inbase'] = True

ensembl_max = ensembl_max.drop(['Gene stable ID version', 'Transcript stable ID', 'Transcript start (bp)', 'Transcript end (bp)', 'length'], axis=1)
ensembl_max.rename(columns={'Gene stable ID':'Gene', 'Transcript stable ID version':'Feature'}, inplace=True)
ensembl_max['blacksign_max'] = True

print(f"Ensembl tables are done [%f sec]" % (time.time() - start_time), end="\n")

for filename in filenames:
	
	print(f"\nStart file {filename} ...", end='\n')
	
	# pickle
	
	start_time = time.time()
	
	with open(f"/dev/datasets/FairWind/_results/dinara/pickles/{filename}.pickle", 'rb') as f:
		main_table = pickle.load(f)
	
	order = main_table.columns.to_list()
	
	print(f"Big data is done [%f sec]" % (time.time() - start_time), end="\n")
	
	# merge
	
	start_time = time.time()
	
	main_table = pd.merge(ensembl_all, main_table, how='right', on=["Gene"])
	main_table = pd.merge(ensembl_max, main_table, how='right', on=["Gene", "Feature"])
	main_table[['blacksign_max', 'blacksign_inbase']] = main_table[['blacksign_max', 'blacksign_inbase']].fillna(False)
	main_table = main_table.loc[~((main_table['blacksign_max'] == False) & (main_table['blacksign_inbase'] == True))]
	main_table = main_table.drop(['blacksign_max', 'blacksign_inbase'], axis=1)
	main_table = main_table[order]
	main_table.sort_values(by=["CHROM","POS"], inplace=True)
	
	print(f"Merging is done [%f sec]" % (time.time() - start_time), end="\n")
	
	# pickling
	
	start_time = time.time()
	    
	with open(f"/dev/datasets/FairWind/_results/dinara/new_pickles/{filename}_light.pickle", 'wb') as f:
		pickle.dump(main_table, f)
	    
	print(f"Pickling is done [%f sec]" % (time.time() - start_time), end="\n")
	
	# write to file
	
	start_time = time.time()
	
	main_table.to_csv(f"/dev/datasets/FairWind/_results/dinara/csv_light/{filename}_light.csv", sep='\t', index=False, mode='w')
	del main_table
	
	print(f"Writing to file is done [%f sec]" % (time.time() - start_time), end="\n")
