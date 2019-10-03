import pandas as pd
 
brother_1 = pd.read_csv("/dev/datasets/FairWind/_results/bowtie/NEW_FILTERED_NAMED/sample-1-3-4_squeezed_named.csv", sep='\t')
brother_2 = pd.read_csv("/dev/datasets/FairWind/_results/bowtie/NEW_FILTERED_NAMED/sample-1-5_squeezed_named.csv", sep='\t')

print(brother_1)
print(brother_2)

pd.merge(brother_1, brother_2, how='inner', on=['CHROM', 'POS', 'REF']).to_csv("/dev/datasets/FairWind/_results/bowtie/merged_brothers.csv", sep='\t', index=False)
