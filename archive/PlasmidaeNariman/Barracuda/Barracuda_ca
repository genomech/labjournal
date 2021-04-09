from Bio import SeqIO
from Bio.Seq import Seq
from multiprocessing import cpu_count, Pool
import sys
import os
import gzip
import subprocess
import re
import pandas
import networkx as nx
import matplotlib.pyplot as plt
import pickle
import json
from pandarallel import pandarallel

#with open("/dev/datasets/FairWind/.cloud/FORGE/data_file.json", "r") as write_file: Connections = json.load(write_file)

#for Barcode, Neighbors in Connections:
	
#print(Connections)
#exit()

pandarallel.initialize(nb_workers=12, verbose=1)

Stats = {"Total": 0, "NotCut": 0, "CutR": 0, "CutF": 0, "CutBoth": 0}

def Bash(Command):

	Shell = subprocess.Popen(Command, shell=True, executable="/bin/bash", stdout=open(os.devnull, 'w'), stderr=subprocess.PIPE)
	Error = Shell.communicate()[1]
	if Shell.returncode != 0: raise OSError(f"Bash fucked up because of architecture proёb")

Input_R1 = "/_home/_temp/nariman_test_R1.fq"
Input_R2 = "/_home/_temp/nariman_test_R2.fq"
Output_R1 = None
Output_R2 = None

SearchLength = 15
BarcodeLengthMismatch = 2
ErrorRate = 0.1
ReadLength = 150

CutFragment = "CATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTGGGCATGCG^NNCGANNACTNNATGNNACGNNCTGNNTCANN$TCGGTACCGAGAACCGGGCAGGTCACGCATCCCCCCCTTCCCTCCCACCCCCTGCCAAGCTCTCCCTCCCAGGATCCTCTCTGGCTCCATCGTAAGCAAACCTTAGAGGTTCTGGCAAGGAGAGAGATGGCTCCAGGAAATGGGGGTGTGTCACCAGATAAGGAATCTGCCTAACAGGAGGTGGGGGTTAGACCCAATATCAGGAGACTAGGAAGGAGGAGGCCTAAGGATGGGGCTTTTCTGTCACCAATCCTGTCCCTAGTGGCCCCACTGTGGGGTGGAGGGGACAGATAAAAGTACCCAGAACCAGAGCCACATTAACCGGCCCTGGGAATATAAGGTGGTCCCAGCTCGGGGACACAGGATCCCTGGAGGCAGCAAACATGCTGTCCTGAAGTGGACATAGGGGCCCGGGTTGGAGGAAGAAGACTAGCTGAGCTCTCGGACCCCTGGAAGATGCCATGACAGGGGGCTGGAAGAGCTAGGGTACCAC^NNTACNNCGANNTAGNNACGNNTTCNNGAGNN$CGCTAGCAACGTAGGAGCGACATTGATTATTGACTAG"

Reverse = "CTAGTCAATAATCAATGTCGCTCCTACGTTGCTAGCG$NNCTCNNGAANNCGTNNCTANNTCGNNGTANN^GTGGTACCCTAGCTCTTCCAGCCCCCTGTCATGGCATCTTCCAGGGGTCCGAGAGCTCAGCTAGTCTTCTTCCTCCAACCCGGGCCCCTATGTCCACTTCAGGACAGCATGTTTGCTGCCTCCAGGGATCCTGTGTCCCCGAGCTGGGACCACCTTATATTCCCAGGGCCGGTTAATGTGGCTCTGGTTCTGGGTACTTTTATCTGTCCCCTCCACCCCACAGTGGGGCCACTAGGGACAGGATTGGTGACAGAAAAGCCCCATCCTTAGGCCTCCTCCTTCCTAGTCTCCTGATATTGGGTCTAACCCCCACCTCCTGTTAGGCAGATTCCTTATCTGGTGACACACCCCCATTTCCTGGAGCCATCTCTCTCCTTGCCAGAACCTCTAAGGTTTGCTTACGATGGAGCCAGAGAGGATCCTGGGAGGGAGAGCTTGGCAGGGGGTGGGAGGGAAGGGGGGGATGCGTGACCTGCCCGGTTCTCGGTACCGA$NNTGANNCAGNNCGTNNCATNNAGTNNTCGNN^CGCATGCCCAGGTGGCACTTTTCGGGGAAATGTGCGCGGAACCCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATG"

CutSeqs = {"1F": CutFragment.split('^')[0][- SearchLength:],
		   "2F": CutFragment.split('$')[1][0: SearchLength],
		   "1R": Reverse.split('$')[0][- SearchLength:],
		   "2R": Reverse.split('^')[1][0: SearchLength],
		   "LF": len(CutFragment.split('^')[1].split('$')[0]),
		   "LR": len(Reverse.split('$')[1].split('^')[0]),
		   "PosF": len(CutFragment.split('^')[0]),
		   "PosR": len(Reverse.split('$')[0])
		   }

print(CutSeqs)

#Command = f"set -o pipefail; cutadapt -j {str(10)} -e {str(ErrorRate)} -g {CutSeqs['1R']} -g {CutSeqs['1F']} -G {CutSeqs['1R']} -G {CutSeqs['1F']} -O {str(SearchLength)} -o /_home/_temp/nariman_test_R1.fq -p /_home/_temp/nariman_test_R2.fq {Input_R1} {Input_R2}"
#Bash(Command)
#exit()

R1 = SeqIO.parse(open(Input_R1, 'r'), "fastq")
R2 = SeqIO.parse(open(Input_R2, 'r'), "fastq")

#R1_write = open(Output_R1, 'w')
#R2_write = open(Output_R2, 'w')

fig = plt.figure(figsize=(12,12))
ax = plt.subplot(111)

G = nx.Graph()

while 1:
	try:
		raise StopIteration
		record_R1 = next(R1)
		record_R2 = next(R2)
		Seq_R1 = record_R1.seq.__str__()
		Seq_R2 = record_R2.seq.__str__()
		Plasmid_R1 = ReadLength - len(Seq_R1)
		Plasmid_R2 = ReadLength - len(Seq_R2)
		Pos_R1 = "F" if ((CutSeqs["PosF"] - BarcodeLengthMismatch) <= Plasmid_R1 <= (CutSeqs["PosF"] + BarcodeLengthMismatch)) else ("R" if ((CutSeqs["PosR"] - BarcodeLengthMismatch) <= Plasmid_R1 <= (CutSeqs["PosR"] + BarcodeLengthMismatch)) else Plasmid_R1)
		Pos_R2 = "F" if ((CutSeqs["PosF"] - BarcodeLengthMismatch) <= Plasmid_R2 <= (CutSeqs["PosF"] + BarcodeLengthMismatch)) else ("R" if ((CutSeqs["PosR"] - BarcodeLengthMismatch) <= Plasmid_R2 <= (CutSeqs["PosR"] + BarcodeLengthMismatch)) else Plasmid_R2)
		Barcode_R1 = "0" if type(Pos_R1) is int else Seq_R1[:32 + BarcodeLengthMismatch]
		Barcode_R2 = "0" if type(Pos_R2) is int else Seq_R2[:32 + BarcodeLengthMismatch]
		
		Stats["Total"] += 1
		
		if (Pos_R1 not in ["F", "R"]) and (Pos_R2 not in ["F", "R"]): Stats["NotCut"] += 1
		elif ((Pos_R1 == "F") and (Pos_R2 not in ["F", "R"])) or ((Pos_R2 == "F") and (Pos_R1 not in ["F", "R"])): Stats["CutF"] += 1
		elif ((Pos_R1 == "R") and (Pos_R2 not in ["F", "R"])) or ((Pos_R2 == "R") and (Pos_R1 not in ["F", "R"])): Stats["CutR"] += 1
		elif (Pos_R1, Pos_R2 == "F", "R") or (Pos_R2, Pos_R1 == "F", "R"): Stats["CutBoth"] += 1
		
		#if Stats["Total"] == 100: raise StopIteration
		
		
		if (Pos_R1 == "F" or Pos_R1 == "R") and (Pos_R2 == "F" or Pos_R2 == "R"):
			Name1 = Pos_R1 + "_" + Barcode_R1
			Name2 = Pos_R2 + "_" + Barcode_R2
			G.add_nodes_from([Name1, Name2])
			if G.has_edge(Name1, Name2): G.edges[Name1, Name2]['weight'] += 1.0
			else: G.add_edge(Name1, Name2, weight=1.0)
		
		continue
		if (Pos_R1 not in ['F', 'R'] and Pos_R2 in ['F', 'R']) or (Pos_R2 not in ['F', 'R'] and Pos_R1 in ['F', 'R']):
			Suffix = " " + str(Pos_R1) + "+" + str(Pos_R2) + "|" + Barcode_R1 + "+" + Barcode_R2
			record_R1.id = record_R1.id + Suffix
			record_R1.description = ""
			record_R2.id = record_R2.id + Suffix
			record_R2.description = ""
			SeqIO.write(record_R1, R1_write, 'fastq')
			SeqIO.write(record_R2, R2_write, 'fastq')
#
	except StopIteration:
		#nx.write_edgelist(G, "test.edgelist.gz")
		#exit()
		#Stats["NotCut"] /= Stats["Total"]
		#Stats["CutF"] /= Stats["Total"]
		#Stats["CutR"] /= Stats["Total"]
		#Stats["CutBoth"] /= Stats["Total"]
		#print(Stats)
		
		#G = nx.read_edgelist("test.edgelist.gz")
		with open('/dev/datasets/FairWind/.cloud/FORGE/test.pickle', 'rb') as f: G = pickle.load(f)
		
		#remove = [edge for edge in G.edges(data=True) if edge[2]["weight"] < 250]
		#G.remove_edges_from(remove)
		#remove = [node for node in G.nodes if not list(G.neighbors(node))]
		
		#G.remove_nodes_from(remove)
		
		#pos = nx.spring_layout(G, k=0.1)
		#nx.draw(G, pos, node_size=2, node_color='red', font_size=8, with_labels=False)

		#plt.tight_layout()
		#plt.savefig("Graph.svg")
		#exit()
		Data = []
		for bc1, bc2 in G.edges:
			Data += [{ "Barcode1": bc1, "Barcode2": bc2"Sum": sum([int(G.edges[Node, x]['weight']) for x in list(G.neighbors(Node))])}]
		Data = pandas.DataFrame(Data).sort_values(by=['Sum'], ascending=False)
		Data.to_csv('/dev/datasets/FairWind/.cloud/FORGE/all_bc.csv')
		
		Data = []
		for Node in list(G.nodes):
			Data += [{ "Barcode": Node, "Sum": sum([int(G.edges[Node, x]['weight']) for x in list(G.neighbors(Node))])}]
		Data = pandas.DataFrame(Data).sort_values(by=['Sum'], ascending=False)
		Data.to_csv('/dev/datasets/FairWind/.cloud/FORGE/all_bc.csv')
		exit()
		Connections = {}
		Target = Data.iloc[:1000,:]
		#Target.to_csv('/dev/datasets/FairWind/.cloud/FORGE/temp.csv')
		#del Target
		#del Data
		
		#Result = pandas.DataFrame(index=list(range(34)), columns=['RA', 'RT', 'RG', 'RC', 'RN', 'FA', 'FT', 'FG', 'FC', 'FN'])
		#Result.fillna(0, inplace=True)
		
		#def Raskidat(Series):
			#Direction = Series["Barcode"][0]
			#Sum = Series["Sum"]
			#Barcode = list(Series["Barcode"][2:])
			#BigList = [[[(Sum if ((Direction == direct) and (item == chars)) else 0) for item in Barcode] for chars in ['A', 'T', 'G', 'C', 'N']] for direct in ['F', 'R']]
			#BigList = [item for sublist in BigList for item in sublist]
			#BigList = [item for sublist in BigList for item in sublist]
			#Header = [[[direct + chars + str(index) for index, item in enumerate(Barcode)] for chars in ['A', 'T', 'G', 'C', 'N']] for direct in ['F', 'R']]
			#Header = [item for sublist in Header for item in sublist]
			#Header = [item for sublist in Header for item in sublist]
			#return pandas.Series(BigList, index=Header)
		
		#Result = []
		#for index, Target in enumerate(pandas.read_csv('/dev/datasets/FairWind/.cloud/FORGE/temp.csv', chunksize=1000)):
			#Target = Target.parallel_apply(Raskidat, axis=1)
			#Target = Target.sum(axis=0)
			#Result += [Target]
			#print(f"Chunk #{index} finished", end='\r')
		
		#Result = pandas.DataFrame(Result)
		#Result = Result.sum(axis=0)
		#Result.to_csv('/dev/datasets/FairWind/.cloud/FORGE/positions.csv')
		##print(Target)
		#exit()
		
		
		#for bc in Target.iterrows():
			#Direction = bc[1]["Barcode"][0]
			#Barcode =  bc[1]["Barcode"][2:]
			#for i in range(34):
				#Result.loc[i,Direction+Barcode[i]] += 1
				#Result.loc[i,Direction+"0"] += 1
		#for direct in ['R', 'F']:
			#for char in ['A', 'T', 'G', 'C', 'N']:
				#Target[direct + char] = Target[direct + char] / Target[direct + "0"]
		
		#Target.to_csv('/dev/datasets/FairWind/.cloud/FORGE/positions.csv')
		#exit()
		
		for bc in Target.iterrows():
			Neighbors = G.neighbors(bc[1]["Barcode"])
			Total = bc[1]["Sum"]
			con = [ { "Neighbor": Node, "Rate": int(G.edges[Node, bc[1]["Barcode"]]['weight']) / Total, "Sum": int(G.edges[Node, bc[1]["Barcode"]]['weight']) } for Node in Neighbors]
			con = pandas.DataFrame(con).sort_values(by=['Rate'], ascending=False)
			con = con.iloc[:10,:]
			neig = {}
			for row in con.iterrows():
				neig[row[1]["Neighbor"]] = [row[1]["Rate"], row[1]["Sum"]]
			Connections[bc[1]["Barcode"]] = neig
		
		#AllSave = [[item] + list(Connections[item]) for item in Connections]
		#AllSave = [item for sublist in AllSave for item in sublist]
		
		with open("/dev/datasets/FairWind/.cloud/FORGE/data_file.json", "w") as write_file: json.dump(Connections, write_file)
		exit()
		#Neighbors = [item for sublist in Neighbors for item in sublist]
		
		#remove = [node for node in G.nodes() if node not in (Target + Neighbors)]
		#G.remove_nodes_from(remove)
		
		#remove = [edge for edge in G.edges(data=True) if edge[2]["weight"] < 200]
		#G.remove_edges_from(remove)
		
		remove = [node for node in G.nodes() if node not in AllSave]
		
		G.remove_nodes_from(remove)
		
		pos = nx.spring_layout(G, k=0.1)
		nx.draw(G, pos, node_size=2, node_color='red', font_size=8, with_labels=False)

		plt.tight_layout()
		plt.savefig("Graph.svg")
		exit()
		
		total = Data['Sum'].sum()
		Data['Rate'] = Data['Sum'].apply(lambda x: x / total)
		#Data['Cumulative rate'] = Data.apply(lambda x: Data['Sum'].iloc[:x.name].sum() / total, axis=1)
		
		edges = sorted(G.edges(data=True), key=lambda t: t[2].get('weight', 1), reverse=True)
		
		edges = pandas.DataFrame(edges, columns=["Barcode1", "Barcode2", "Weight"])
		edges["Weight"] = edges["Weight"].apply(lambda x: int(x["weight"]))
		edges.to_csv("/dev/datasets/FairWind/_temp/BarcodesEdges.tsv", sep='\t', index=False)
		Data.to_csv("/dev/datasets/FairWind/_temp/BarcodesTable.tsv", sep='\t', index=False)
		
		Freq = Data.groupby('Sum').count()
		Freq = Freq.reset_index()
		Freq.groupby(["Sum", pandas.cut(Freq["Sum"], 300, right=False)])['Barcode'].sum().unstack(fill_value=0).sum() 

		
		print("Basta, diablo!")
		break
