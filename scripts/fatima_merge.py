from lib.blister import *

ncbi = pd.read_csv("/dev/datasets/FairWind/_results/Fatima/gene_id_names.csv", sep='\t') #, header=None, names=['gene_id', 'ncbi_id']) 

table = pd.read_csv("/dev/datasets/FairWind/_results/Fatima/cuffdiff/gene_exp.diff", sep='\t')
ncbi.rename(columns={'Cuff ID' : 'gene_id'}, inplace=True)
table = pd.merge(ncbi, table, how='right', on=['gene_id'])
table.sort_values(by=['log2(fold_change)'], ascending=False, inplace=True)
table['NCBI ID'].fillna(value='-', inplace=True)

table = table[(table['q_value'] < 0.05)] # & ((table['log2(fold_change)'] < -1) & (table['log2(fold_change)'] > float('-inf')))]

table.sort_values(by=['log2(fold_change)'], axis=0, ascending=True, inplace=True)
table.rename(columns={'value_1':'WT', 'value_2':'Mut'}, inplace=True)
table.drop(columns=['gene', 'test_id', 'locus', 'sample_1', 'sample_2', 'status', 'test_stat', 'p_value', 'significant'], inplace=True)
table.fillna(value='-', inplace=True)

print(table)

table.to_csv("/dev/datasets/FairWind/_results/Fatima/sayeeda_AllGenesDiff_named.csv", sep=',', index=False)
