from lib.blister import *

ncbi = pd.read_csv("/dev/datasets/FairWind/_results/Andre/gene_id.csv", sep='\t', header=None, names=['gene_id', 'ncbi_id']) 

table = pd.read_csv("/dev/datasets/FairWind/_results/Andre/cuffdiff/gene_exp.diff", sep='\t')

table = pd.merge(ncbi, table, how='right')
table.sort_values(by=['log2(fold_change)'], ascending=False, inplace=True)
table['ncbi_id'].fillna(value='-', inplace=True)

table = table[(table['q_value'] < 0.05) & ((table['log2(fold_change)'] < -1) & (table['log2(fold_change)'] > float('-inf')))]

table.sort_values(by=['log2(fold_change)'], axis=0, ascending=True, inplace=True)
table.rename(columns={'value_1':'Day', 'value_2':'Night'}, inplace=True)
table.drop(columns=['gene', 'test_id', 'locus', 'sample_1', 'sample_2', 'status', 'test_stat', 'p_value', 'significant'], inplace=True)

print(Blister.GitHubTable(table))

#table.to_csv("/dev/datasets/FairWind/_results/Andre/genes_plus_diff.csv", sep=',', index=False)
