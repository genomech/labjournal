from lib.cuffcake import *

table = pd.read_csv("/dev/datasets/FairWind/_results/bowtie/1-9_NF1_vep.csv", sep='\t')
table['chrom'] = table["Location"].apply(lambda x: x.split(':')[0])
table['begin'] = table["Location"].apply(lambda x: x.split(':')[1].split('-')[0])
table['end'] = table["Location"].apply(lambda x: x.split(':')[1].split('-')[0])

reference = pd.read_csv("/dev/datasets/FairWind/_db/rn6.bed", sep='\t', header=None, names=['chrom', 'begin', 'end', 1, 2, 3, 4, 5, 6, 7, 8, 9])

Cuffcake.LociIntersection(table, reference, table_format='three', reference_format='three')
