from lib.blister import *
from _Liebe_banished import *

Blister.Logo("Test bleat!")

with Blister.Timestamp("Parsing Table") as start_time:
	table = pd.read_csv("/dev/datasets/FairWind/_results/Andre/cuffdiff_4gr_v2/gene_exp.diff", sep='\t')

with Blister.Timestamp("Parsing Reference") as start_time:
	reference = pd.read_csv("/dev/datasets/FairWind/_db/norm_names.csv", sep='\t')
	reference.rename(columns={'txStart' : 'begin', 'txEnd' : 'end'}, inplace=True)
	reference = reference[['name', 'name2', 'chrom', 'begin', 'end']]

with Blister.Timestamp("Expression") as start_time:
	exp = Cuffcake.ExpressionPerSample("/dev/datasets/FairWind/_results/Andre/cuffdiff_4gr_v2", genes_id=None)

with Blister.Timestamp("Prepare Table") as start_time:
	table["significant"] = table["significant"].apply(lambda x: True if x == 'yes' else False)
	for col in ['value_1', 'value_2', 'log2(fold_change)', 'test_stat', 'p_value', 'q_value']:
		table[col] = pd.to_numeric(table[col], downcast='float', errors='raise')
	table = table[(table['q_value'] < 0.05) & (~((table['sample_1'] == 'SD-NIGHT') & (table['sample_2'] == 'L27-DAY'))) & (~((table['sample_1'] == 'SD-DAY') & (table['sample_2'] == 'L27-NIGHT')))]

with Blister.Timestamp("Loci Intersection") as start_time:
	new_table = Cuffcake.LociIntersection(table, reference, table_format='one', reference_format='three')

new_table = new_table[['gene_id', 'name_ref', 'name2_ref', 'sample_1', 'sample_2', 'status', 'value_1', 'value_2', 'log2(fold_change)', 'q_value',  'chrom', 'start', 'end', 'start_ref', 'end_ref', 'intersection_ref']]
new_table = pd.merge(exp, new_table, how='right', right_on='gene_id', left_index=True)
new_table = new_table[['gene_id', 'name_ref', 'name2_ref', 'sample_1', 'sample_2', 'status', 'value_1', 'value_2', 'log2(fold_change)', 'q_value', 'chrom', 'start', 'end', 'start_ref', 'end_ref', 'intersection_ref', 'MB_AFR_001_S1_sorted.bam, SD-DAY, 0', 'MB_AFR_002_S2_sorted.bam, SD-DAY, 1', 'MB_AFR_003_S3_sorted.bam, SD-DAY, 2', 'MB_AFR_004_S4_sorted.bam, SD-NIGHT, 0', 'MB_AFR_005_S5_sorted.bam, SD-NIGHT, 1', 'MB_AFR_006_S6_sorted.bam, SD-NIGHT, 2', 'MB_AFR_007_S7_sorted.bam, L27-DAY, 0', 'MB_AFR_008_S8_sorted.bam, L27-DAY, 1', 'MB_AFR_009_S9_sorted.bam, L27-DAY, 2', 'MB_AFR_010_S10_sorted.bam, L27-NIGHT, 0', 'MB_AFR_011_S11_sorted.bam, L27-NIGHT, 1', 'MB_AFR_012_S12_sorted.bam, L27-NIGHT, 2']]
new_table.fillna('-', inplace=True)
new_table.sort_values(by=['sample_1', 'sample_2', 'log2(fold_change)', 'value_1', 'value_2'], ascending=False, inplace=True)
new_table.to_csv("/dev/datasets/FairWind/_results/Andre/4gr_data_v3.csv", sep='\t', index=False, float_format="%.5f")

