from lib.blister import *
import numpy as np

def the_thread(block, output_dir):
	columns_list = ['AF', 'AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'AA_AF', 'EA_AF', 'gnomAD_AF', 'gnomAD_AFR_AF', 'gnomAD_AMR_AF', 'gnomAD_ASJ_AF', 'gnomAD_EAS_AF', 'gnomAD_FIN_AF', 'gnomAD_NFE_AF', 'gnomAD_OTH_AF', 'gnomAD_SAS_AF']
	index, input_filename = block
	Blister.Sleep(max_level=70.0, interval=1, index=index)
	table = pd.read_csv(input_filename, sep='\t')
	for col in columns_list:
		table[col] = pd.to_numeric(table[col], errors='raise', downcast='float')
	table['TOTAL_AF'] = table[columns_list].apply(np.nanmax, axis=1)
	table = table[((table['IMPACT'] == 'HIGH') | (table['IMPACT'] == 'MODERATE')) & ((table['TOTAL_AF'] < 0.01) | table['TOTAL_AF'].isna())]
	print(table)
	
results = Blister.EachFile("Tomsk Filter Tables", ["/dev/datasets/FairWind/_results/bowtie/tables_severe_filtered/sample-1-1_Filtered.csv"], "/dev/datasets/FairWind/_results/bowtie/NEW_FILTERED", THREADS_NUM = 6)(the_thread)()
