import sys
import pandas as pd
import numpy as np
import string
import re
from Bio import SeqIO
from lib.blister import *

C_INPUT_FILES = ["/dev/datasets/FairWind/_results/November/Illuminaless/*1_Illuminaless.fq.gz"]
C_DIR_NAME = "/dev/datasets/FairWind/_results/November/November_BGG_100k"
C_GENOME = '-genome-'
C_SEQUENCES = np.array([
	['-illumina-', 'AGATCGGAAG'],
	['-bridge--gatc-', 'GCTGAGGGATC'],
	['-bridge-', 'GCTGAGG'],
	['-gatc--egdirb-', 'GATCCCTCAGC'],
	['-egdirb-', 'CCTCAGC'],
	['-blunt--gac-', 'CAGTGGCGAC'],
	['-blunt-', 'CAGTGGC'],
	['-gtc--tnulb-', 'GTCGCCACTG'],
	['-tnulb-', 'GCCACTG']
	])
C_KMER_SIZE = 12
C_MAX = 100000

def thread0(filename, output_dir):

	seqs = pd.DataFrame(C_SEQUENCES, columns=['name', 'seq'])
	genome = C_GENOME
	output_file0 = Blister.Output(filename[1], output_dir, "BGG-100k", "fastq", rewrite=True, index=filename[0])
	main_list = []
	main_table = pd.DataFrame(columns=['count', 'mask'])
	total = 0

	with Blister.Read(filename[1], 'rt', index=filename[0]) as input_file, open(output_file0, "w") as output_handle:
		with Blister.Timestamp(f"PARSING", index=filename[0]) as start_time:
			for record in SeqIO.parse(input_file, "fastq"):
				line = record.seq.__str__()
				for index, row in seqs.iterrows(): line = line.replace(row['seq'], row['name'])
				for it in range(1, C_KMER_SIZE): line = re.sub("(^|-)[ATGCN]{%s}($|-)" % (it), "--[%s]--" % it, line)
				line = re.sub(r"[ATGCN]+", genome, line)
				if line == "-bridge--gatc--genome-":
					SeqIO.write(record, output_handle, "fastq")
					total += 1
				if total == C_MAX: break

Blister.Logo("Analysis of Sequences v3.1 [blister]")

input_filenames = Blister.Input(C_INPUT_FILES)
if not input_filenames: exit()
output_dir = Blister.Dir(C_DIR_NAME, create=True)
if not output_dir: exit()

with Blister.Threading(f"ANALYSIS") as pool:
	pool.map(functools.partial(thread0, output_dir=output_dir), enumerate(input_filenames))

