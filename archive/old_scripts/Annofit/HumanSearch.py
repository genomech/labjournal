#!/bin/python3

# TODO Z-score ExAC

PLI_DOM = 0.9

Legend = {'PVS1': 0,
		  'PS1': 2, 'PS2': 3, 'PS3': 4, 'PS4': 5,
		  'PM1': 10, 'PM2': 11, 'PM3': 12, 'PM4': 13, 'PM5': 14, 'PM6': 15,
		  'PP1': 16, 'PP2': 17, 'PP3': 18, 'PP4': 19, 'PP5': 20,
		  'BA1': 1,
		  'BS1': 6, 'BS2': 7, 'BS3': 8, 'BS4': 9,
		  'BP1': 21, 'BP2': 22, 'BP3': 23, 'BP4': 24, 'BP5': 25, 'BP6': 26, 'BP7': 27}

def sortir(element): return Legend[element]

Patogenic = ['PVS1',
			 'PS1', 'PS2', 'PS3', 'PS4',
			 'PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6',
			 'PP1', 'PP2', 'PP3', 'PP4', 'PP5']

import pandas as pd
import sys

input_filename = sys.argv[1]
output_filename = sys.argv[2]

output_format = output_filename.split('.')[-1]
assert (output_format in ['xlsx', 'csv']), 'Output format must be .xlsx or .csv'

header = ['Comment', 'RANK', 'HGMD', 'avsnp150', 'Chr', 'Start', 'End', 'GeneName_comp', 'Ref', 'Alt', 'VCF_QUAL', 'Func_comp', 'ExonicFunc_comp', 'PopFreqMax', 'pc_H', 'pc_M', 'pc_L', 'pc_U', 'Significance_comp', 'VCF_GT', 'pLI', 'OMIM_Phenotypes', 'OMIM_linked', 'OMIM_dominance', 'Compound', 'regsnp_splicing_site', 'CLNDN']
exonic_func = ['stoploss', 'stopgain', 'startloss', 'frameshift insertion', 'frameshift deletion']
significance = ['Likely pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic', 'Pathogenic', 'association', '_association', 'risk_factor', '_Affects', '_risk_factor']
ncrna = ['ncRNA_splicing', 'ncRNA_intronic', 'ncRNA_exonic']

data = pd.read_csv(input_filename, sep='\t', na_values='.')
data['Comment'] = ""
data['RANK'] = data.apply(lambda x: [], axis=1)

##  CODES

# Get Existing Codes
def getCodes(line):
	if line['PVS1'] != line['PVS1']: return line['RANK']
	for item in list(Legend.keys()):
		if int(line[item]) == 1: line['RANK'] += [item]
	return line['RANK']

data['RANK'] = data[['RANK'] + list(Legend)].apply(getCodes, axis=1)

# Population: set BS1 (AD 0.01%, AR 0.5%, XD 0.03%, XR 0.5%) and PM2 (AD 0.01%, AR 0.5%, XD 0.01%, XR 0.3%)
# pLI > 0.9 == D
def func_BS1_PM2(line):
	if (line['pLI'] != line['pLI']) and (line['OMIM_dominance'] != line['OMIM_dominance']): return line['RANK']
	
	line['pLI'] = max([float(x) for x in str(line['pLI']).split(';')])
	if (line['pLI'] >= PLI_DOM) or (line['OMIM_dominance'] == 'D') or (line['OMIM_dominance'] == 'RD'):
		if ((line['OMIM_linked'] == 'X') or (line['Chr'] == 'chrX')):
			if (line['PopFreqMax'] >= 0.03): return line['RANK'] + ['BS1']
			if (0 < line['PopFreqMax'] < 0.01): return line['RANK'] + ['PM2']
		if (line['OMIM_linked'] == 'A') or ((line['Chr'] != 'chrX') and (line['Chr'] != 'chrY')):
			if (line['PopFreqMax'] >= 0.01): return line['RANK'] + ['BS1']
			if (0 < line['PopFreqMax'] < 0.01): return line['RANK'] + ['PM2']
	
	if (line['pLI'] < PLI_DOM) or (line['OMIM_dominance'] == 'R'):
		if ((line['OMIM_linked'] == 'X') or (line['Chr'] == 'chrX')):
			if (line['PopFreqMax'] >= 0.5): return line['RANK'] + ['BS1']
			if (0 < line['PopFreqMax'] < 0.3): return line['RANK'] + ['PM2']
		if (line['OMIM_linked'] == 'A') or ((line['Chr'] != 'chrX') and (line['Chr'] != 'chrY')):
			if (line['PopFreqMax'] >= 0.5): return line['RANK'] + ['BS1']
			if (0 < line['PopFreqMax'] < 0.5): return line['RANK'] + ['PM2']
	return line['RANK']

data['RANK'] = data[['RANK', 'OMIM_linked', 'OMIM_dominance', 'pLI', 'Chr', 'PopFreqMax']].apply(func_BS1_PM2, axis=1)

# Impact: set BP4 (>=3 Low Prediction Score) & PP3 (>=3 High Prediction Score)
data['RANK'] = data[['RANK', 'pc_L', 'pc_H']].apply(lambda x: x['RANK'] + ['BP4'] if (x['pc_L'] >= 3) else (x['RANK'] + ['PP3'] if (x['pc_H'] >= 3) else x['RANK']), axis=1)

# Function: set PVS1 (stoploss, stopgain, startloss, frameshift, splicing)
data['RANK'] = data[['RANK', 'ExonicFunc_comp', 'Func_comp']].apply(lambda x: x['RANK'] + ['PVS1'] if ( any([(y == 'splicing') for y in str(x['Func_comp']).replace(';', ',').split(',')]) or any([(y in exonic_func) for y in str(x['ExonicFunc_comp']).split(',')])) else x['RANK'], axis=1)

# RANK Prepare

def Ranker(item):
	item = [x[:-1] for x in list(set(item))]
	cPVS, cPS, cPM, cPP, cBA, cBS, cBP = item.count('PVS'), item.count('PS'), item.count('PM'), item.count('PP'), item.count('BA'), item.count('BS'), item.count('BP')
	
	patogenic_1 = (cPVS >= 1) and ((cPS >= 1) or (cPM >= 2) or ((cPM >= 1) and (cPP >= 1)) or (cPP >= 2))
	patogenic_2 = (cPS >= 2)
	patogenic_3 = (cPS >= 1) and ((cPM >= 3) or ((cPM >= 2) and (cPP >= 2)) or ((cPM >= 1) and (cPP >= 4)))
	patogenic_s = patogenic_1 or patogenic_2 or patogenic_3
	
	likely_patogenic_1 = (cPVS >= 1) and (cPM >= 1)
	likely_patogenic_2 = (cPS >= 1) and (cPM >= 1)
	likely_patogenic_3 = (cPS >= 1) and (cPP >= 2)
	likely_patogenic_4 = (cPM >= 3)
	likely_patogenic_5 = (cPM >= 2) and (cPP >= 2)
	likely_patogenic_6 = (cPM >= 1) and (cPP >= 4)
	likely_patogenic_s = likely_patogenic_1 or likely_patogenic_2 or likely_patogenic_3 or likely_patogenic_4 or likely_patogenic_5 or likely_patogenic_6
	
	benign_1 = (cBA >= 1)
	benign_2 = (cBS >= 2)
	benign_s = benign_1 or benign_2
	
	likely_benign_1 = (cBS >= 1) and (cBP >= 1)
	likely_benign_2 = (cBP >= 2)
	likely_benign_s = likely_benign_1 or likely_benign_2
	
	uncertain_1 = not (patogenic_s or likely_patogenic_s or likely_benign_s or benign_s)
	uncertain_2 = (patogenic_s or likely_patogenic_s) and (benign_s or likely_benign_s)
	uncertain_s = uncertain_1 or uncertain_2
	
	if uncertain_s: return 'Uncertain'
	if patogenic_s: return 'Pathogenic'
	if benign_s: return 'Benign'
	if likely_patogenic_s: return 'Likely pathogenic'
	if likely_benign_s: return 'Likely benign'
	return '.'

data['RANK_summary'] = data['RANK'].apply(Ranker)
data['RANK'] = data['RANK'].apply(lambda x: ','.join(sorted(list(set(x)), key=sortir)))

data['Significance_comp'] = data[['Significance_comp', 'RANK_summary']].apply(lambda x: str(x['Significance_comp']) + ',' + str(x['RANK_summary']), axis=1)
data['RANK'] = data['Significance_comp'].apply(lambda x: sum([(y in significance) for y in str(x).split(',')]))
data.drop_duplicates(["Chr", "Start", "End", "Ref", "Alt"], inplace=True)

## FILTERS
omim_filter = data['OMIM_Phenotypes'].apply(lambda x: type(x) is not float)
hgmd_filter = data['HGMD'].apply(lambda x: type(x) is not float)
population_filter = data['PopFreqMax'] < 0.03 # remove BA1 (<= 3%)
impact_filter = (data['pc_H'] >= 3)
significance_filter = data['Significance_comp'].apply(lambda x: False if (type(x) is not str) else any([(y in significance) for y in x.split(',')]))
splicing_filter = data['Func_comp'].apply(lambda x: False if (type(x) is not str) else any([(y == 'splicing') for y in x.replace(';', ',').split(',')]))
ncrna_filter = data['Func_comp'].apply(lambda x: False if (type(x) is not str) else any([(y in ncrna) for y in x.replace(';', ',').split(',')]))
exonicfunc_filter = data['ExonicFunc_comp'].apply(lambda x: False if (type(x) is not str) else any([(y in exonic_func) for y in x.split(',')]))

# Filters Apply
data = data[population_filter & (hgmd_filter | impact_filter | significance_filter | exonicfunc_filter | splicing_filter | (ncrna_filter & omim_filter))]

# Compound, zygosity & dominance
pli_filter = data['pLI'].apply(lambda x: False if x != x else any([float(y) >= PLI_DOM for y in str(x).split(';')]))
dominance_filter = data['OMIM_dominance'].apply(lambda x: False if x != x else ((x == 'RD') or (x == 'D')))
homozygosity_filter = data['VCF_GT'].apply(lambda x: False if type(x) is not float else x == 'HOMO')
noinfo_filter = data['pLI'].apply(lambda x: x != x) & data['OMIM_dominance'].apply(lambda x: x != x)
def compound_count(gene): return data['GeneName_comp'].apply(lambda x: gene in x.replace(';', ',').split(',')).value_counts()[True]
data['Compound'] = data['GeneName_comp'].apply(lambda x: ','.join([str(compound_count(gene)) for gene in x.replace(';', ',').split(',')]))
compound_filter = data['Compound'].apply(lambda x: any([int(y) > 1 for y in x.split(',')]))
data = data[compound_filter | pli_filter | dominance_filter | homozygosity_filter | noinfo_filter]

# Save
data = data[header]
data.fillna('.', inplace=True)
if output_format == 'xlsx': data.to_excel(output_filename, index=False)
elif output_format == 'csv': data.to_csv(output_filename, sep='\t', index=False)
