import pysam
from lib.blister import *
import copy

def the_thread(block, output_dir):
	
	first = True
	dick = [{}, {}, {}]
	counter = [0, 0, 0]
	index, input_filename = block
	output_filename = Blister.Output(input_filename, output_dir, "", "bam", rewrite=True, index=index)
	output_fasta = Blister.Output(input_filename, output_dir, "", "fa", rewrite=True, index=index)
	samfile = pysam.AlignmentFile(input_filename, "r")
	pairedreads = pysam.AlignmentFile(output_filename, "wb", template=samfile)
	buf = None
	for read in samfile:
		target = 0 if (0 < read.reference_start < 400) else (1 if (3000 < read.reference_start < 3500) else 2)
		if first: buf = (read.seq, read, ((70, 145) if (not read.is_reverse) else (-145, -70)))
		else:
			coords = ((70, 145) if (not read.is_reverse) else (-145, -70))
			try:
				dick[target][read.seq[coords[0]:coords[1]]][1] += 1
			except KeyError:
				dick[target][read.seq[coords[0]:coords[1]]] = [read, 1, {} ]
			try:
				dick[target][buf[0][buf[2][0]:buf[2][1]]][1] += 1
			except KeyError:
				dick[target][buf[0][buf[2][0]:buf[2][1]]] = [buf[1], 1, {} ]
			
			try:
				dick[target][read.seq[coords[0]:coords[1]]][2][buf[0]] += 1
			except KeyError:
				dick[target][read.seq[coords[0]:coords[1]]][2][buf[0]] = 1
			try:
				dick[target][buf[0][buf[2][0]:buf[2][1]]][2][read.seq] += 1
			except KeyError:
				dick[target][buf[0][buf[2][0]:buf[2][1]]][2][read.seq] = 1
			
		counter[target] += 1
		first = not first
	fasta = open(output_fasta, 'wt')
	for t in range(3):
		cock = {}
		for key, value in dick[t].items():
			value[0].query_name = value[0].query_name + "__(" + str(value[1]) + ")"
			try:
				cock[value[1]][0] += 1
			except KeyError:
				cock[value[1]] = [1, value[0], value[2]]
		s1 = pd.Series(cock).sort_index()
		s2 = copy.deepcopy(s1)
		for index2, value2 in s2.items():
			sus = 0
			for index1, value1 in s1.items():
				if index1 >= index2: sus += index1 * value1[0]
			s2[index2] = sus * 100 / counter[t]
		
		s1.name = "num"
		s2.name = "cumulative"
		q = s1.to_frame()
		q["cumulative"] = s2
		q["cumulative"] = pd.to_numeric(q["cumulative"], errors='coerce')
		q = q.nsmallest(20, "cumulative") # top 30?
		
		for item in q["num"].items():
			rock = {}
			for key, value in item[1][2].items():
				try:
					rock[value][0] += 1
				except KeyError:
					rock[value] = [1, key]
			s10 = pd.Series(rock).sort_index()
			s20 = copy.deepcopy(s10)
			for index2, value2 in s20.items():
				sus = 0
				counter0 = 0
				for index1, value1 in s10.items(): counter0 += index1 * value1[0]
				for index1, value1 in s10.items():
					if index1 >= index2: sus += index1 * value1[0]
				s20[index2] = sus * 100 / counter0
			s10.name = "num"
			s20.name = "cumulative"
			q0 = s10.to_frame()
			q0["cumulative"] = s20
			q0 = q0[q0["cumulative"] < 95.0]
			q["num"][item[0]][2] = q0["num"].to_list()
		for item in q["num"].items(): 
			pairedreads.write(item[1][1])
			coords1 = ((70, 145) if (not item[1][1].is_reverse) else (-145, -70))
			seq = item[1][1].seq[coords1[0]:coords1[1]]
			fasta.write(">" + item[1][1].query_name + '\n')
			fasta.write(seq + '\n')
			#q["num"][item[0]][1] = q["num"][item[0]][1].seq
		#q.to_csv(output_table, sep='\t')
	fasta.close()

Blister.EachFile("Clusterization", ["/dev/datasets/FairWind/Pavel/sorted_cut/a*.bam"], "/dev/datasets/FairWind/Pavel/clusterized/", THREADS_NUM = 4)(the_thread)()
