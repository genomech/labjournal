from lib.blister import *

coverage = pd.read_csv("/dev/datasets/FairWind/_results/Andre/cuffdiff/genes.fpkm_tracking", sep='\t')
coverage.drop(columns=['tracking_id', 'class_code', 'nearest_ref_id'], inplace=True, axis=1)
print(coverage)
