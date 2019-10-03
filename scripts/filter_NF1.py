from lib.blister import *
import numpy as np

def the_thread(block, output_dir):
	
	columns_list = ['AF', 'AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'AA_AF', 'EA_AF', 'gnomAD_AF', 'gnomAD_AFR_AF', 'gnomAD_AMR_AF', 'gnomAD_ASJ_AF', 'gnomAD_EAS_AF', 'gnomAD_FIN_AF', 'gnomAD_NFE_AF', 'gnomAD_OTH_AF', 'gnomAD_SAS_AF']
	
	
	
	index, input_filename = block
	output_filename = Blister.Output(input_filename, output_dir, "squeezed", "csv", rewrite=True, index=index)
	if not output_filename: return

	Blister.Sleep(max_level=70.0, interval=1, index=index)

	with Blister.Timestamp("Processing file", input_filename, output_filename, index=index) as start_time:
		table = pd.read_csv(input_filename, sep='\t', low_memory=False)
		for col in columns_list:
			table[col] = pd.to_numeric(table[col], errors='coerce', downcast='float')
		table['TOTAL_AF'] = table[columns_list].apply(np.nanmax, axis=1)
		table['TOTAL_AF'] = table['TOTAL_AF'].apply(lambda x: x if x == x else 0)
		#table = table[((table['IMPACT'] == 'HIGH') | (table['IMPACT'] == 'MODERATE')) & ((table['TOTAL_AF'] < 0.01) | table['TOTAL_AF'].isna())]
		
		#table = pd.merge(genemap, table, how='right', on=['SYMBOL'])
		
		#table = table[order]
		order = list(table.columns.values)
		table = table.groupby(table['Location'])
	
		new_table = pd.DataFrame()
		for col in order:
			new_table[col] = table[col].apply(set).apply(list).apply(lambda x: [t for t in x if str(t) != 'nan']).apply(lambda x: (x[0] if ((len(x) == 1) and (type(x) == type(list()))) else (float('nan') if ((not x) and (type(x) == type(list()))) else ', '.join([str(t) for t in x]))))
		del table
		new_table.sort_values(by=['TOTAL_AF'], axis=0, ascending=True, inplace=True)
		new_table.to_csv(output_filename, sep='\t', index=False)
		del new_table

results = Blister.EachFile("Tomsk Filter Tables", ["/dev/datasets/FairWind/_results/bowtie/1-9_NF1_vep.csv"], "/dev/datasets/FairWind/_results/bowtie/")(the_thread)()
