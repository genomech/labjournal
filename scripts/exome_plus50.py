from lib.blister import *
import copy

Blister.Logo("Exome +50!")

contour = 50

exome = pd.read_csv("/dev/datasets/FairWind/_db/MedExome_hg19_capture_targets.sorted.bed", sep='\t', header=None, names=['chrom', 'start', 'end', 'trash'])
for col in ['start', 'end']: exome[col] = pd.to_numeric(exome[col], downcast='unsigned', errors='raise')
exome.drop(columns=['trash'], inplace=True)

genome = pd.read_csv("/dev/datasets/FairWind/_db/hg19/hg19.bed", sep='\t', header=None, names=['chrom', 'start', 'end'])
for col in ['start', 'end']: genome[col] = pd.to_numeric(genome[col], downcast='unsigned', errors='raise')
genome.set_index(['chrom'], drop=True, append=False, inplace=True, verify_integrity=True)

new_exome = copy.deepcopy(exome)
new_exome["start"] = new_exome["start"].apply(lambda x: x - contour)
new_exome["end"] = new_exome["end"].apply(lambda x: x + contour)

def theta(x):
	if x['start'] < 0: x['start'] = 0
	end = genome.at[x['chrom'], 'end']
	if x['end'] > end: x['end'] = end
	return x

new_exome = new_exome.apply(theta, axis=1)
new_exome.to_csv("/dev/datasets/FairWind/_results/CNTN6/exome_plus50.bed", sep='\t', index=False, header=False)
