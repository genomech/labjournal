from lib.blister import *
import copy

Blister.Logo("Not Exome!")

exome = pd.read_csv("/dev/datasets/FairWind/_db/MedExome_hg19_capture_targets.sorted.bed", sep='\t', header=None, names=['chrom', 'start', 'end', 'trash'])
for col in ['start', 'end']: exome[col] = pd.to_numeric(exome[col], downcast='unsigned', errors='raise')
exome.drop(columns=['trash'], inplace=True)

genome = pd.read_csv("/dev/datasets/FairWind/_db/hg19/hg19.bed", sep='\t', header=None, names=['chrom', 'start', 'end'])
for col in ['start', 'end']: genome[col] = pd.to_numeric(genome[col], downcast='unsigned', errors='raise')
genome.set_index(['chrom'], drop=True, append=False, inplace=True, verify_integrity=True)

exome_before = copy.deepcopy(exome)
exome_before.index = exome_before.index.map(lambda x: x + 1)
exome_before.columns = exome_before.columns.map(lambda x: x + "_before")
exome_before.drop(columns=['start_before'], inplace=True)

exome_after = copy.deepcopy(exome)
exome_after.index = exome_after.index.map(lambda x: x - 1)
exome_after.columns = exome_after.columns.map(lambda x: x + "_after")
exome_after.drop(columns=['end_after'], inplace=True)

exome = pd.merge(exome_before, exome, how='right', left_index=True, right_index=True)
exome = pd.merge(exome_after, exome, how='right', left_index=True, right_index=True)

exome.fillna('-1', inplace=True)
for col in ['start_after', 'end_before']: exome[col] = exome[col].apply(int)

def theta(x):
	if x['chrom_before'] != x['chrom']: x['end_before'] = 0
	if x['chrom_after'] != x['chrom']: x['start_after'] = genome[genome.index == x['chrom']]['end'][0]
	return x

exome = exome.apply(theta, axis=1)

results = []

for row in exome.iterrows():
	if row[1]['chrom_before'] != row[1]['chrom']: results += [ [ row[1]['chrom'], row[1]['end_before'], row[1]['start'] ] ]
	results += [ [ row[1]['chrom'], row[1]['end'], row[1]['start_after'] ] ]

not_exome = pd.DataFrame(results, columns=['chrom', 'start', 'end'])
not_exome['trash'] = '.'
not_exome.to_csv("/dev/datasets/FairWind/_db/NOT_MedExome_hg19.bed", sep='\t', index=False, header=False)
