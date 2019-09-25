from lib.blister import *
import intervals as I

class Node: 
	def __init__(self, key): 
		self.left = None
		self.right = None
		self.val = key

def insert(root, node):
	if root is None: root = node 
	else:
		if root.val['interval'].lower < node.val['interval'].lower: 
			if root.right is None: root.right = node 
			else: insert(root.right, node) 
		else: 
			if root.left is None: root.left = node 
			else: insert(root.left, node)

def search(root, value): 
	if root is None: return None
	else:
		if root.val['interval'].is_connected(value): return root.val['ncbi_id']
		else:
			if root.val['interval'].lower < value.lower: return search(root.right, value) 
			else: return search(root.left, value)

def inorder(root): 
	if root: 
		inorder(root.left) 
		print(root.val)
		inorder(root.right) 

Blister.Logo("Sayeeda's Gene Finder")

with Blister.Timestamp("LOAD TRANSCRIPTOME") as start_time:
	transcriptome = pd.read_csv("/dev/datasets/FairWind/_results/Andre/cuffmerge/merged.gtf", sep='\t', header=None, names=['chrom', 'hz1', 'hz2', 'begin', 'end', 'hz3', 'hz4', 'transcript_id', 'gene_id'])
	transcriptome.drop(columns=['hz1', 'hz2', 'hz3', 'hz4'], inplace=True, axis=0)
	transcriptome['transcript_id'] = transcriptome['gene_id'].apply(lambda x: x.split(" ")[3][1:-2])
	transcriptome['gene_id'] = transcriptome['gene_id'].apply(lambda x: x.split(" ")[1][1:-2])
	for col in ['begin', 'end']:
		transcriptome[col] = pd.to_numeric(transcriptome[col], downcast='integer', errors='raise')
	
	length = transcriptome.shape[0]
	intervals = []
	for it in transcriptome.iterrows():
		intervals += [I.IntInterval.closed(it[1]['begin'], it[1]['end'])]
		assert not I.IntInterval.closed(it[1]['begin'], it[1]['end']).empty
		Blister.ProgressBar(it[0] / length, start_time)
		
	transcriptome['interval'] = intervals
	transcriptome.drop(columns=['begin', 'end'], inplace=True, axis=0)
	del intervals
	Blister.Erase()

with Blister.Timestamp("LOAD GENOME") as start_time:
	genome = pd.read_csv("/dev/datasets/FairWind/_db/rn6.bed", sep='\t', header=None, names=['chrom', 'begin', 'end', 'ncbi_id', 'hz1', 'hz2', 'hz3', 'hz4', 'hz5', 'hz6', 'hz7', 'hz8'])
	genome.drop(columns=['hz1', 'hz2', 'hz3', 'hz4', 'hz5', 'hz6', 'hz7', 'hz8'], inplace=True, axis=0)
	for col in ['begin', 'end']:
		genome[col] = pd.to_numeric(genome[col], downcast='integer', errors='raise')

	chroms = list(set(genome['chrom'].to_list()))
	sorted_genome = dict()
	done = 0
	for chrom in chroms:
		table = genome[genome['chrom'] == chrom]
		trigger = False
		sorted_genome[chrom] = None
		for it in table.iterrows():
			it[1]['interval'] = I.IntInterval.closed(it[1]['begin'], it[1]['end'])
			it[1].drop(labels=['begin', 'end', 'chrom'], inplace=True)
			if trigger: insert(sorted_genome[chrom], Node(it[1]))
			else:
				sorted_genome[chrom] = Node(it[1])
				trigger = True
		del table
		done += 1
		Blister.ProgressBar(done / len(chroms), start_time)
	Blister.Erase()

with Blister.Timestamp("FIND GENES") as start_time:
	length = transcriptome.shape[0]
	ncbi = []
	for trans in transcriptome.iterrows():
		if trans[1]['chrom'] in sorted_genome.keys():
			ncbi += [search(sorted_genome[trans[1]['chrom']], trans[1]['interval'])]
			Blister.ProgressBar(trans[0] / length, start_time)
		else:
			ncbi += [None]
	transcriptome['ncbi_id'] = ncbi
	Blister.Erase()

with Blister.Timestamp("PROCESS TABLE & SAVE") as start_time:
	transcriptome = transcriptome['ncbi_id'].groupby(transcriptome['gene_id'])
	transcriptome = transcriptome.apply(set).apply(list).apply(lambda x: [i for i in x if i]).apply(lambda x: x[0] if x else None)
	transcriptome.to_csv("/dev/datasets/FairWind/_results/Andre/gene_id.csv", sep='\t', index=True)
