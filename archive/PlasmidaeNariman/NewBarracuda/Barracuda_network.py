from Bio import SeqIO
from Bio.Seq import Seq
from multiprocessing import cpu_count, Pool
import sys
import os
import gzip
import re
import pandas
import networkx as nx
import pickle
import json
from pandarallel import pandarallel
from copy import deepcopy as dc

Input_R1 = "/dev/datasets/ngs_data/Battulin_BGI_20200918/demult/AS_A23_Cut_1.fq.gz"
Input_R2 = "/dev/datasets/ngs_data/Battulin_BGI_20200918/demult/AS_A23_Cut_2.fq.gz"
EdgeList = "/dev/datasets/ngs_data/Battulin_BGI_20200918/demult/Edgelist.list"
BigDataFile = "/dev/datasets/ngs_data/Battulin_BGI_20200918/demult/BigData.tsv"
StatsFile = "/dev/datasets/ngs_data/Battulin_BGI_20200918/demult/Stats.tsv"

CutSeqs = {'1F': 'CCACCTGGGCATGCG', '2F': 'TCGGTACCGAGAACC', 'LF': 32, 'PosF': 91, '1R': 'CCTACGTTGCTAGCG', '2R': 'GTGGTACCCTAGCTC', 'LR': 32, 'PosR': 37}

Stats = {"Absolute": {"Total": 0, "No cut": 0, "Forward cut": 0, "Reverse cut": 0, "FR cut": 0, "Strange (FF | RR)": 0 } }

ReadLength = 150
BarcodeLengthMismatch = 2
Threads = 10

pandarallel.initialize(nb_workers=Threads, verbose=1)

Stream = {
	"R1": SeqIO.parse(gzip.open(Input_R1, 'rt'), "fastq"),
	"R2": SeqIO.parse(gzip.open(Input_R2, 'rt'), "fastq"),
	"BigData": open(BigDataFile, 'wt')
	}

G = nx.Graph()

while True:
	try:
		
		# Prepare
		Data = {"R1": {}, "R2": {}}
		for Mate in ["R1", "R2"]:
			Data[Mate]["Record"] = next(Stream[Mate])
			Data[Mate]["Sequence"] = dc(Data[Mate]["Record"].seq.__str__())
			Data[Mate]["Length"] = ReadLength - len(Data[Mate]["Sequence"])
			Data[Mate]["Position"] = "F" if ((CutSeqs["PosF"] - BarcodeLengthMismatch) <= Data[Mate]["Length"] <= (CutSeqs["PosF"] + BarcodeLengthMismatch)) else ("R" if ((CutSeqs["PosR"] - BarcodeLengthMismatch) <= Data[Mate]["Length"] <= (CutSeqs["PosR"] + BarcodeLengthMismatch)) else None)
			Data[Mate]["Barcode"] = None if Data[Mate]["Position"] is None else Data[Mate]["Sequence"][:CutSeqs[f"L{Data[Mate]['Position']}"]+BarcodeLengthMismatch]
			Data[Mate]["Record"] = None
		Stream["BigData"].write('\t'.join([str(Data[i][j]) for i in Data.keys() for j in Data[i].keys()]) + '\n')
		
		PositionCheck = lambda Data, Param, Value1, Value2: ((Data["R1"][Param] == Value1) and (Data["R2"][Param] == Value2)) or ((Data["R2"][Param] == Value1) and (Data["R1"][Param] == Value2))
		
		# Graph
		if (Data["R1"]["Position"] is not None) and (Data["R2"]["Position"] is not None):
			Names = {"R1": f"{Data['R1']['Position']}_{Data['R1']['Barcode']}", "R2": f"{Data['R2']['Position']}_{Data['R2']['Barcode']}"}
			G.add_nodes_from([Names["R1"], Names["R2"]])
			if G.has_edge(Names["R1"], Names["R2"]): G.edges[Names["R1"], Names["R2"]]['weight'] += 1.0
			else: G.add_edge(Names["R1"], Names["R2"], weight=1.0)
		
		# Stats
		Stats["Absolute"]["Total"] += 1
		Stats["Absolute"]["No cut"] += 1 if PositionCheck(Data, "Position", None, None) else 0
		Stats["Absolute"]["Forward cut"] += 1 if PositionCheck(Data, "Position", 'F', None) else 0
		Stats["Absolute"]["Reverse cut"] += 1 if PositionCheck(Data, "Position", 'R', None) else 0
		Stats["Absolute"]["FR cut"] += 1 if PositionCheck(Data, "Position", 'F', 'R') else 0
		Stats["Absolute"]["Strange (FF | RR)"] += 1 if (PositionCheck(Data, "Position", 'F', 'F') or PositionCheck(Data, "Position", 'R', 'R')) else 0
		print(Stats["Absolute"]["Total"], end='\r')
		Data = None
		#raise StopIteration
	except StopIteration:
		Stream["BigData"].close()
		break

# Stats

Stats["Rates, %"] = {key: value / Stats["Absolute"]["Total"] * 100 for key, value in Stats["Absolute"].items()}
Stats = pandas.DataFrame(Stats)
Stats.to_csv(StatsFile, sep='\t')

print("Do edgelist...", end='\n')
pickle.dump(G, open(f"{EdgeList}.pickle", 'wb'))
nx.write_edgelist(G, open(EdgeList, 'wb'))
