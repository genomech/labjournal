import pandas as pd
import gzip
import re
import sys

data = pd.read_csv("/dev/datasets/FairWind/_db/MedExome_hg19_capture_targets.sorted.bed", sep='\t', header=None)
data[data[0].apply(lambda x: re.fullmatch("^chr((\d+)|X|Y)$", str(x)) is not None)][[0, 1, 2]].to_csv("/dev/datasets/FairWind/_db/MedExome_hg19_capture_targets.pipeline.bed", index=False, sep='\t', header=None)

exit()

chunksize = 10 ** 7
input_filename, output_filename = sys.argv[1], sys.argv[2]

with gzip.open(input_filename, 'rt') as input_file, open(output_filename, 'wt') as output_file:
	
	x = True
	
	while x: 
		line = input_file.readline()
		if re.match("^##contig\=<ID\=.*,length\=\d+,assembly\=.*>$", line) is not None:
			if re.match("^##contig\=<ID\=((\d+)|X|Y),length\=\d+,assembly\=.*>$", line) is not None: line = line.replace("<ID=", "<ID=chr")
			else: continue
		output_file.write(line)
		x = line[0:2] == '##'
	
	for chunk in pd.read_csv(input_file, chunksize=chunksize, sep='\t', header=None):
		chunk = chunk[chunk[0].apply(lambda x: re.fullmatch("^((\d+)|X|Y)$", str(x)) is not None)]
		chunk[0] = chunk[0].apply(lambda x: "chr" + str(x))
		output_file.write(chunk.to_csv(index=False, sep='\t', header=None)) 
