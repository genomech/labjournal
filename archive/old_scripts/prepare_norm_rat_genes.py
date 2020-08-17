from lib.blister import *

table = pd.read_csv("/dev/datasets/FairWind/_db/norm_names.csv", sep='\t')
table.drop(columns=['#bin', 'chrom', 'strand', 'txStart', 'txEnd', 'exonEnds', 'score', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'cdsStartStat', 'cdsEndStat', 'exonFrames'], inplace=True)
table['name'] = table['name'].apply(lambda x: x[0:-2])
table.rename(columns={'name' : 'NCBI ID', 'name2' : 'Gene Symbol'}, inplace=True)

old_table = pd.read_csv("/dev/datasets/FairWind/_results/Fatima/gene_id.csv", header=None, names=['Cuff ID', 'NCBI ID'], sep='\t')
table = pd.merge(old_table, table, how='right', on=['NCBI ID'])
table = table[table['Cuff ID'] == table['Cuff ID']]

table.to_csv("/dev/datasets/FairWind/_results/Fatima/gene_id_names.csv", sep='\t', index=False)
