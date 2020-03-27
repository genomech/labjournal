import pandas as pd

header = ['HGMD', 'avsnp150', 'Chr', 'Start', 'End', 'GeneName_comp', 'Ref', 'Alt', 'VCF_QUAL', 'Func_comp', 'ExonicFunc_comp', 'PopFreqMax', 'pc_H', 'pc_M', 'pc_L', 'pc_U', 'Significance_comp', 'VCF_GT', 'pLI', 'OMIM_Phenotypes', 'OMIM_linked', 'OMIM_dominance', 'regsnp_splicing_site', 'CLNDN']

files = ['TAF3_S1_L001_recalibrated_daemonically_annotated.csv', 'TAF4_S2_L002_recalibrated_daemonically_annotated.csv']

for f in files:
	data = pd.read_csv("/dev/datasets/FairWind/_results/CNTN6/tables/" + f, sep='\t', na_values='.')
	data = data[header]
	data = data[((data['PopFreqMax'] < 0.05) & (data['pc_H'] > 1)]
	
	#data = data[data['OMIM_Phenotypes'].apply(lambda x: type(x) is str)]
	#keywords = ['bile', 'chol', 'renal', 'reno', 'kidney', 'mucoviscido'] #  'syndrome', 'disease'
	#data["keywords"] = data['OMIM_Phenotypes'].apply(lambda x: any(ele in x.lower() for ele in keywords) if (type(x) is str) else False)
	#data = data[data["keywords"]]

	data.fillna('.', inplace=True)
	data.to_excel('/dev/datasets/FairWind/_results/CNTN6/filtered_2/' + f.split('.')[0] + ".xlsx", index=False)
