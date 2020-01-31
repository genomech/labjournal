import pandas as pd
import json
import numpy as np

with open("otherinfo.json", 'r') as k: other_info = json.load(k)
with open("order.list", 'r') as f: order = f.read().splitlines()
merge_func = lambda x: ','.join(list(set(x.to_list())))
prediction_merge = lambda x: ''.join(x.to_list())
data = pd.read_csv("/dev/datasets/FairWind/_results/case/128_S5_cover_10-200_primitive.vcf.hg19_multianno.txt", sep='\t')

# hgmd + omim
hgmd_data = pd.read_csv("hgmd.csv", sep='\t')
omim_data = pd.read_csv("omim.csv", sep='\t')
data = pd.merge(hgmd_data, data, how='right', on=["Chr", "Start"])
data = pd.merge(omim_data, data, how='right', on=["Gene.refGene"])
data.sort_values(by=["Chr", "Start"], inplace=True)

# other
data.rename(columns=other_info, inplace=True)
data["Func_comp"] = data[["Func.knownGene", "Func.refGene", "Func.ensGene"]].apply(merge_func, axis=1)
data["ExonicFunc_comp"] = data[["ExonicFunc.knownGene", "ExonicFunc.refGene", "ExonicFunc.ensGene"]].apply(merge_func, axis=1)
data["GeneName_comp"] = data[["Gene.knownGene", "Gene.refGene"]].apply(merge_func, axis=1)
data["Prediction_comp"] = data[["SIFT_pred", "LRT_pred", "MutationTaster_pred", "MutationAssessor_pred", "FATHMM_pred", "PROVEAN_pred", "MetaSVM_pred", "MetaLR_pred", "M-CAP_pred", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred", "LR_pred", "regsnp_disease"]].apply(prediction_merge, axis=1)
data["Significance_comp"] = data[["CLNSIG", "InterVar_automated"]].apply(merge_func, axis=1)
data.fillna(".", inplace=True)
data["Filter_1"] = np.nan
data["Filter_2"] = np.nan
data = data[order]
data.to_csv("/dev/datasets/FairWind/_results/case/128_S5_daemonically_annotated.csv", sep='\t', index=False)

