import pysam
import pandas as pd

primers = pd.DataFrame([
	['L2F', 'AGCACGTGGGGTGAGAG', 75, 91, 'F'],
	['L2R', 'TGCTCAAGTAGACCTAATGTGG', 279, 300, 'R'],
	['R2F', 'ACCTCATGATCCAAGGGTACCTCC', 4207, 4230, 'F'],
	['R2R', 'TCGCCTGGAATCCTCCAGCT', 4416, 4435, 'R'],
	['M2F', 'AGAGTGTGGCTGGGTACCTG', 3123, 3142, 'F'],
	['M2R', 'CTTGGCCACACAGGTGTAGTT', 3329, 3349, 'R']
	], columns=['name', 'seq', 'start', 'end', 'orient'])

reverse_barcodes = {'GCAT' : 'ATGC', 'CATG' : 'CATG', 'ATGC' : 'GCAT', 'TGCA' : 'TGCA'}

# in the table below columns barcodes are reversal.

barcodes = pd.DataFrame([
	['ATGC', False, 1, 2, 3],
	['CATG', 4, 5, 6, 7],
	['GCAT', 8, 9, 10, 11],
	['TGCA', 12, 13, 14, 15],
	['CGTA', 16, 17, 18, 19]
	], columns=['forward', 'GCAT', 'CATG', 'ATGC', 'TGCA'])

samfile = pysam.AlignmentFile("head.sam", "r")

table_cols = ['name', 'read1', 'read2', 'primer1', 'primer2', 'barcode1', 'barcode2', 'animal']

main_table = pd.DataFrame(columns=table_cols)

total = 0
first = True
_buffer = None

for read in samfile:
	print(read.query_sequence)
	break
	if first:
		_buffer = pd.Series(index=table_cols)
		_buffer['name'] = read.query_name
		_buffer['read1'] = read
		primer1 = []
		for s in primers.iterrows():
			if s[1]['seq']
	first = not first
	total += 1
	if total == 10: break

samfile.close()
