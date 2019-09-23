from lib.blister import *

ncbi = pd.read_csv("/dev/datasets/FairWind/_results/Fatima/gene_id.csv", sep='\t', header=None, names=['gene_id', 'ncbi_id']) 

table = pd.read_csv("/dev/datasets/FairWind/_results/Fatima/genes_diff.csv", sep='\t')

table = pd.merge(ncbi, table, how='right')
table.sort_values(by=['log2_fold_change'], ascending=False, inplace=True)
table['ncbi_id'].fillna(value='-', inplace=True)
print(table)

table.to_csv("/dev/datasets/FairWind/_results/Fatima/genes_diff_ncbi.csv", sep='\t', index=False)
