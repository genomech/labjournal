from lib.blister import *

Blister.Logo("Andre's Beasts")

tracking = pd.read_csv("/dev/datasets/FairWind/_results/Andre/cuffdiff/genes.read_group_tracking", sep='\t')
samples = pd.read_csv("/dev/datasets/FairWind/_results/Andre/cuffdiff/read_groups.info", sep='\t')
samples.rename(columns={'replicate_num' : 'replicate'}, inplace=True)
samples.drop(columns=['total_mass', 'norm_mass', 'internal_scale', 'external_scale'], inplace=True, axis=1)
tracking = pd.merge(samples, tracking, on=['replicate', 'condition'], how='right')

lst = ['XLOC_002266', 'XLOC_011074', 'XLOC_017099', 'XLOC_002641', 'XLOC_017100', 'XLOC_005619', 'XLOC_017752', 'XLOC_000883', 'XLOC_006483', 'XLOC_000169', 'XLOC_000668', 'XLOC_013338', 'XLOC_017751', 'XLOC_012607', 'XLOC_015082', 'XLOC_017361', 'XLOC_006079', 'XLOC_016089', 'XLOC_015909', 'XLOC_016036', 'XLOC_007781', 'XLOC_017757', 'XLOC_015224', 'XLOC_014627', 'XLOC_001103', 'XLOC_001784', 'XLOC_006086', 'XLOC_009955']

tracking = tracking[tracking['tracking_id'].isin(lst)]
tracking.drop(columns=['raw_frags', 'internal_scaled_frags', 'external_scaled_frags', 'effective_length', 'status'], inplace=True, axis=1)
tracking.sort_values(by=['tracking_id', 'condition', 'replicate'], ascending=True, inplace=True)
cols = tracking.columns
tracking = tracking.groupby('tracking_id').apply(pd.Series.to_list).apply(pd.DataFrame, columns=cols).apply(pd.DataFrame.drop, columns=['condition', 'replicate', 'tracking_id'], axis=1).to_dict()

table = 1

for it in tracking.keys():
	tracking[it] = tracking[it].rename(columns={'FPKM' : it})
	if type(table) == type(pd.DataFrame()): table = pd.merge(tracking[it], table, on=['file'], how='right')
	else: 
		table = tracking[it]

table['file'] = table['file'].apply(lambda x: x[7:10])

table = table.transpose()
table.columns = ['01 [D]', '02 [D]', '03 [D]', '07 [D]', '09 [D]', '04 [N]', '05 [N]', '06 [N]', '08 [N]', '10 [N]', '11 [N]', '12 [N]']
table.drop(labels='file', axis=0, inplace=True)
table.to_csv("/dev/datasets/FairWind/.cloud/core/labjournal/labjournal/scripts_results/andre_s_beasts.csv", sep=',', index=True)
