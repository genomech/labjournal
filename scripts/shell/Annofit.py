# perl $annovar_dir/table_annovar.pl $input_vcf $annovar_dir/humandb -buildver $genome_assembly -protocol knownGene,ensGene,refGene,abraom,AFR.sites.2015_08,ALL.sites.2015_08,AMR.sites.2015_08,ASN.sites.2012_04,avgwas_20150121,avsift,avsnp150,cadd13,cg69,clinvar_20190305,cosmic70,dann,dbnsfp35c,dbscsnv11,EAS.sites.2015_08,eigen,esp6500_all,EUR.sites.2015_08,exac03,fathmm,gene4denovo201907,gerp++,gme,gnomad211_genome,gwava,hrcr1,icgc21,intervar_20180118,kaviar_20150923,ljb26_all,mcap13,mitimpact24,MT_ensGene,nci60,popfreq_all_20150413,regsnpintron,revel,SAS.sites.2015_08,snp142 --operation g,g,g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f --remove --vcfinput --thread $threads;

__version__ = "0.02"
__author__ = "regnveig"

import collections
import numpy as np
import os
import pandas as pd
import sys

current_dir = os.path.dirname(sys.argv[0])
input_filename = sys.argv[1]
output_filename = sys.argv[2]

# func
def df_squeeze(dataframe):
	squeezed = pd.Series()
	for col in dataframe.columns.to_list():
		squeezed[col] = ';'.join(sorted([str(x) for x in dataframe[col].to_list() if str(x) != 'nan']))
		if squeezed[col] == '': squeezed[col] = float('NaN')
	return squeezed

merge_func = lambda x: ','.join(list(set(x.to_list())))
prediction_merge = lambda x: ''.join(x.to_list())

# const
PC_header = ["pc_H", "pc_M", "pc_L", "pc_U"]
vcf1_header = ["VCF_GT", "VCF_DP", "VCF_AD", "VCF_RO", "VCF_QR", "VCF_AO", "VCF_QA", "VCF_GL"]
prediction_set = ["SIFT_pred", "LRT_pred", "MutationTaster_pred", "MutationAssessor_pred", "FATHMM_pred", "PROVEAN_pred", "MetaSVM_pred", "MetaLR_pred", "M-CAP_pred", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred", "LR_pred", "regsnp_disease"]
numbered_cols = ["PopFreqMax", "AF_male", "AF_female", "non_neuro_AF_popmax", "non_cancer_AF_popmax", "VCF_DP"]
other_info = {"Otherinfo1": "trash_1", "Otherinfo2": "VCF_QUAL", "Otherinfo3": "trash_2", "Otherinfo4": "VCF_Chr", "Otherinfo5": "VCF_Start", "Otherinfo6": "trash_3", "Otherinfo7": "VCF_Ref", "Otherinfo8": "VCF_Alt", "Otherinfo9": "trash_4", "Otherinfo10": "trash_5", "Otherinfo11": "VCF_Info", "Otherinfo12": "VCF_DataHeader", "Otherinfo13": "VCF_DataValues"}

# load data
with open(current_dir + "/annofit_data/order.list", 'r') as f: order = f.read().splitlines()
hgmd_data = pd.read_csv(current_dir + "/annofit_data/hgmd.csv", sep='\t')
omim_data = pd.read_csv(current_dir + "/annofit_data/omim_pli.csv", sep='\t')
data = pd.read_csv(input_filename, sep='\t')

# hgmd + omim
data = pd.merge(hgmd_data, data, how='right', on=["Chr", "Start"])
data["Gene.refGene"] = data["Gene.refGene"].apply(lambda x: ';'.join(sorted(list(set(str.split(x, ';'))))))
genez = [str.split(x, ';') for x in data["Gene.refGene"][data["Gene.refGene"].apply(lambda x: x.find(';') != -1)].unique()]
for item in genez: omim_data = omim_data.append(df_squeeze(omim_data[omim_data["Gene.refGene"].apply(lambda x: x in item)]), ignore_index=True)
omim_data.dropna(axis=0, how='all', inplace=True)
data = pd.merge(omim_data, data, how='right', on=["Gene.refGene"])
data.rename(columns=other_info, inplace=True)

# vcf data
gotcha = data["VCF_DataValues"].apply(lambda x: str.split(str(x), ':'))
gotcha = gotcha.apply(lambda x: pd.Series(x, index=vcf1_header) if len(x) == len(vcf1_header) else pd.Series(x + ([float('NaN')] * 7), index=vcf1_header))
gotcha["VCF_GT"] = gotcha["VCF_GT"].apply(lambda x: x.replace('|', '/')).apply(lambda x: "-" if (x == "0/0") else ("HOMO" if (len(set(x.split('/'))) == 1) else x))
gotcha["VCF_DP"].fillna('.', inplace=True)
data = pd.concat([data, gotcha], axis=1)

# other
data["Func_comp"] = data[["Func.knownGene", "Func.refGene", "Func.ensGene"]].apply(merge_func, axis=1)
data["ExonicFunc_comp"] = data[["ExonicFunc.knownGene", "ExonicFunc.refGene", "ExonicFunc.ensGene"]].apply(merge_func, axis=1)
data["GeneName_comp"] = data[["Gene.knownGene", "Gene.refGene"]].apply(merge_func, axis=1)
data["Significance_comp"] = data[["CLNSIG", "InterVar_automated"]].apply(merge_func, axis=1)

data["Prediction_comp"] = data[prediction_set].apply(prediction_merge, axis=1)
predict = data["Prediction_comp"].apply(lambda x: pd.Series([x.count('D') + x.count('H') + x.count('A'), x.count('P') + x.count('M'), x.count('T') + x.count('B') + x.count('N') + x.count('L'), x.count('.') + x.count('U')], index=PC_header))
data = pd.concat([data, predict], axis=1)

for col in numbered_cols: data[col] = data[col].apply(lambda x: -1 if x == '.' else x)
data.fillna(".", inplace=True)

data["Filter_1"] = np.nan
data["Filter_2"] = np.nan
data["Filter_3"] = np.nan
data["Filter_4"] = np.nan
data = data[order]

# sort
data.sort_values(by=["Chr", "Start"], inplace=True)

data.to_csv(output_filename, sep='\t', index=False)
