#!/bin/python3

__version__ = "0.02"
__author__ = "Zanthia"

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
		squeezed[col] = ';'.join(sorted([x for x in dataframe[col].to_list() if str(x) != 'nan']))
		if squeezed[col] == '': squeezed[col] = float('NaN')
	return squeezed

merge_func = lambda x: ','.join(list(set(x.to_list())))
prediction_merge = lambda x: ''.join(x.to_list())

# const
vcf1_header = ["VCF_GT", "VCF_DP", "VCF_AD", "VCF_RO", "VCF_QR", "VCF_AO", "VCF_QA", "VCF_GL"]
prediction_set = ["SIFT_pred", "LRT_pred", "MutationTaster_pred", "MutationAssessor_pred", "FATHMM_pred", "PROVEAN_pred", "MetaSVM_pred", "MetaLR_pred", "M-CAP_pred", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred", "LR_pred", "regsnp_disease"]
numbered_cols = ["PopFreqMax", "AF_male", "AF_female", "non_neuro_AF_popmax", "non_cancer_AF_popmax", "VCF_DP"]
other_info = {"Otherinfo1": "trash_1", "Otherinfo2": "VCF_QUAL", "Otherinfo3": "trash_2", "Otherinfo4": "VCF_Chr", "Otherinfo5": "VCF_Start", "Otherinfo6": "trash_3", "Otherinfo7": "VCF_Ref", "Otherinfo8": "VCF_Alt", "Otherinfo9": "trash_4", "Otherinfo10": "trash_5", "Otherinfo11": "VCF_Info", "Otherinfo12": "VCF_DataHeader", "Otherinfo13": "VCF_DataValues"}

# load data
with open(current_dir + "/annofit_data/order.list", 'r') as f: order = f.read().splitlines()
hgmd_data = pd.read_csv(current_dir + "/annofit_data/hgmd.csv", sep='\t')
omim_data = pd.read_csv(current_dir + "/annofit_data/omim_splinted.csv", sep='\t')
data = pd.read_csv(input_filename, sep='\t')

# hgmd + omim
data = pd.merge(hgmd_data, data, how='right', on=["Chr", "Start"])
data["Gene.refGene"] = data["Gene.refGene"].apply(lambda x: ';'.join(sorted(str.split(x, ';'))))
genez = [str.split(x, ';') for x in data["Gene.refGene"][data["Gene.refGene"].apply(lambda x: x.find(';') != -1)].unique()]
for item in genez: omim_data = omim_data.append(df_squeeze(omim_data[omim_data["Gene.refGene"].apply(lambda x: x in item)]), ignore_index=True)
omim_data.dropna(axis=0, how='all', inplace=True)
data = pd.merge(omim_data, data, how='right', on=["Gene.refGene"])
data.rename(columns=other_info, inplace=True)

# vcf data
gotcha = data["VCF_DataValues"].apply(lambda x: str.split(str(x), ':'))
gotcha = gotcha.apply(lambda x: pd.Series(x, index=vcf1_header) if len(x) == len(vcf1_header) else pd.Series(x + ([float('NaN')] * 7), index=vcf1_header))
gotcha["VCF_GT"] = gotcha["VCF_GT"].apply(lambda x: x.replace('|', '/'))
gotcha["VCF_DP"].fillna('.', inplace=True)
data = pd.concat([data, gotcha], axis=1)

# sort
data.sort_values(by=["Chr", "Start"], inplace=True)

# other
data["Func_comp"] = data[["Func.knownGene", "Func.refGene", "Func.ensGene"]].apply(merge_func, axis=1)
data["ExonicFunc_comp"] = data[["ExonicFunc.knownGene", "ExonicFunc.refGene", "ExonicFunc.ensGene"]].apply(merge_func, axis=1)
data["GeneName_comp"] = data[["Gene.knownGene", "Gene.refGene"]].apply(merge_func, axis=1)
data["Prediction_comp"] = data[prediction_set].apply(prediction_merge, axis=1)
data["Significance_comp"] = data[["CLNSIG", "InterVar_automated"]].apply(merge_func, axis=1)

for col in numbered_cols: data[col] = data[col].apply(lambda x: -1 if x == '.' else x)
data.fillna(".", inplace=True)

data["Filter_1"] = np.nan
data["Filter_2"] = np.nan
data["Filter_3"] = np.nan
data["Filter_4"] = np.nan
data = data[order]

data.to_csv(output_filename, sep='\t', index=False)
