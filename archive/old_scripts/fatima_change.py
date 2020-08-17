from lib.blister import *
import matplotlib.pyplot as plt
import numpy as np

ncbi = pd.read_csv("/dev/datasets/FairWind/_results/Fatima/gene_id_names.csv", sep='\t')

genes_diff = pd.read_csv("/dev/datasets/FairWind/_results/Andre/cuffdiff_4gr/gene_exp.diff", sep='\t')
genes_diff["significant"] = genes_diff["significant"].apply(lambda x: True if x == 'yes' else False)
for col in ['value_1', 'value_2', 'log2(fold_change)', 'test_stat', 'p_value', 'q_value']:
	genes_diff[col] = pd.to_numeric(genes_diff[col], downcast='float', errors='raise')

genes_diff = genes_diff[(genes_diff['q_value'] < 0.05) & (~((genes_diff['sample_1'] == 'SD-NIGHT') & (genes_diff['sample_2'] == 'L27-DAY'))) & (~((genes_diff['sample_1'] == 'SD-DAY') & (genes_diff['sample_2'] == 'L27-NIGHT')))]

ncbi.rename(columns={'Cuff ID' : 'gene_id'}, inplace=True)
genes_diff = pd.merge(ncbi, genes_diff, how='right', on=['gene_id'])
genes_diff.fillna(value='-', inplace=True)

genes_diff.sort_values(by=['sample_1', 'sample_2', 'log2(fold_change)', 'value_1', 'value_2'], axis=0, ascending=False, inplace=True)
genes_diff.drop(columns=['status', 'test_id', 'gene', 'test_stat', 'p_value', 'significant'], inplace=True)

#genes_diff.to_csv("/dev/datasets/FairWind/_results/Andre/diff_4gr.csv", sep=",", index=False)

print(set(genes_diff['Gene Symbol'].to_list()))
exit()

genes_diff = pd.read_csv("file:///dev/datasets/FairWind/_results/Andre/cuffdiff_4gr/gene_exp.diff", sep='\t')
genes_diff["sample_1"] = 'wt'
genes_diff["sample_2"] = 'mut'
genes_diff["significant"] = genes_diff["significant"].apply(lambda x: True if x == 'yes' else False)
for col in ['value_1', 'value_2', 'log2_fold_change', 'test_stat', 'p_value', 'q_value']:
	genes_diff[col] = pd.to_numeric(genes_diff[col], downcast='float', errors='raise')

genes_diff = genes_diff[(genes_diff['q_value'] < 0.05) & ((genes_diff['log2_fold_change'] < -1) & (genes_diff['log2_fold_change'] > float('-inf')))]
#genes_diff = genes_diff[(genes_diff['q_value'] < 0.05) & (genes_diff['log2_fold_change'] == float('-inf'))]

#genes_diff["log2_fold_change"] = np.absolute(genes_diff["log2_fold_change"])
genes_diff.sort_values(by=['log2_fold_change'], axis=0, ascending=True, inplace=True)
genes_diff.rename(columns={'value_1':'wt', 'value_2':'mut'}, inplace=True)
genes_diff.drop(columns=['sample_1', 'sample_2', 'status', 'test_stat', 'p_value', 'significant'], inplace=True)
genes_diff.to_csv("/dev/datasets/FairWind/_results/Fatima/genes_minus_diff.csv", sep="\t", index=False)
print(genes_diff)
