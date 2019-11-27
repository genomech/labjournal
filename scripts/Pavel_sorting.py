import pysam
import pandas as pd
from lib.blister import *

C_SHORT = 8
C_MAX = 0
C_POS_MISMATCH = 4

primers = pd.DataFrame([
	['L2', 'AGCACGTGGGGTGAGAG', 75, 91, 'F'],
	['L2', 'TGCTCAAGTAGACCTAATGTGG', 279, 300, 'R'],
	['R2', 'ACCTCATGATCCAAGGGTACCTCC', 4207, 4230, 'F'],
	['R2', 'TCGCCTGGAATCCTCCAGCT', 4416, 4435, 'R'],
	['M2', 'AGAGTGTGGCTGGGTACCTG', 3123, 3142, 'F'],
	['M2', 'CTTGGCCACACAGGTGTAGTT', 3329, 3349, 'R']
	], columns=['name', 'seq', 'start', 'end', 'orient'])

primers['shorten'] = primers.apply(lambda x: x[1][:C_SHORT], axis=1)

reverse_barcodes = {'GCAT' : 'ATGC', 'CATG' : 'CATG', 'ATGC' : 'GCAT', 'TGCA' : 'TGCA'}

# in the table below columns barcodes are reversal.

barcodes = pd.DataFrame([
	[False, 1, 2, 3],
	[4, 5, 6, 7],
	[8, 9, 10, 11],
	[12, 13, 14, 15],
	[16, 17, 18, 19]
	], index=['ATGC', 'CATG', 'GCAT', 'TGCA', 'CGTA'], columns=['GCAT', 'CATG', 'ATGC', 'TGCA'])

# ------ CODE ------

Blister.Logo("Pavel Animal Sorting")

def processing(block):
	AllData = {}
	AllData['result'] = {}
	AllData['animal'] = 0
	for key in block.keys():
		for primer in primers.iterrows():
			data = {}
			pos_left = block[key].query_sequence.find(primer[1]['seq'])
			if pos_left != -1:
				pos_right = len(block[key].query_sequence) - (pos_left + len(primer[1]['seq']))
				data['read_pos'] = pos_left if pos_left < pos_right else pos_right
				data['primer_name'] = primer[1]['name']
				data['orient'] = 'F' if pos_left < pos_right else 'R'
				data['pos'] = (block[key].pos + 1) if pos_left < pos_right else (block[key].pos + block[key].query_length - 6)
				data['right_orient'] = data['orient'] == primer[1]['orient']
				data['dpos'] = (primer[1]['start'] if pos_left < pos_right else primer[1]['end']) - data['pos']
				data['barcode'] = block[key].query_sequence[(pos_left - 4):pos_left] if pos_left < pos_right else block[key].query_sequence[(- pos_right):(- pos_right + 4)]
				break
			data = None
		AllData['result'][key] = data
	
	if type(AllData['result'][1]) == type(dict()) and type(AllData['result'][2]) == type(dict()):
		r1 = AllData['result'][1] if AllData['result'][1]['pos'] < AllData['result'][2]['pos'] else AllData['result'][2]
		r2 = AllData['result'][1] if AllData['result'][1]['pos'] > AllData['result'][2]['pos'] else AllData['result'][2]
		if (r1['barcode'] in barcodes.index) and (r2['barcode'] in barcodes.columns) and (r1['primer_name'] == r2['primer_name']) and (abs(r1['dpos']) < C_POS_MISMATCH) and (abs(r2['dpos']) < C_POS_MISMATCH) and r1['right_orient'] and r2['right_orient']:
			AllData['animal'] = barcodes[r2['barcode']][r1['barcode']]
			for key in block.keys(): block[key].qname = block[key].qname + "__" + AllData['result'][key]['primer_name'] + AllData['result'][key]['orient'] + ":" + AllData['result'][key]['barcode']
			return (AllData['animal'], block)
	return (0, block)

results = {}
for i in range(0,20): results[i] = []

total = 0
first = True
_buffer = {1 : None, 2 : None }

samfile = pysam.AlignmentFile("/dev/datasets/FairWind/Pavel/seq.sam", "r")

with Blister.Timestamp("SORT") as start_time:
	for key in results.keys(): results[key] = pysam.AlignmentFile(f"/dev/datasets/FairWind/Pavel/sorted/animal_{key}.sam", "w", template=samfile)
	for read in samfile:
		if first: _buffer[1] = read
		else:
			_buffer[2] = read
			assert _buffer[1].query_name == _buffer[2].query_name
			result = processing(_buffer)
			results[result[0]].write(result[1][1])
			results[result[0]].write(result[1][2])
			Blister.ProgressBar(total / 51044527, start_time)
			total += 1
		first = not first
		if C_MAX != 0:
			if total == C_MAX: break
	Blister.Erase()
	for key in results.keys(): result[key].close()
	samfile.close()

print(f"Total pairs: {total}")
