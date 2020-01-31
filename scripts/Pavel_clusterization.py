import pysam
from lib.blister import *
import copy

def the_thread(block, output_dir):
	dick = {}
	cock = {}
	counter = 0
	index, input_filename = block
	output_filename = Blister.Output(input_filename, output_dir, "", "bam", rewrite=True, index=index)
	output_table = Blister.Output(input_filename, output_dir, "table", "txt", rewrite=True, index=index)
	samfile = pysam.AlignmentFile(input_filename, "r")
	pairedreads = pysam.AlignmentFile(output_filename, "wb", template=samfile)
	
	for read in samfile:
		try:
			dick[read.seq][1] += 1
		except KeyError:
			dick[read.seq] = [read, 1]
		counter += 1
	
	for key, value in dick.items():
		value[0].query_name = value[0].query_name + "__(" + str(value[1]) + ")"
		pairedreads.write(value[0])
		try:
			cock[value[1]] += 1
		except KeyError:
			cock[value[1]] = 1
	s1 = pd.Series(cock).sort_index()
	s2 = copy.deepcopy(s1)
	for index2, value2 in s2.items():
		sus = 0
		for index1, value1 in s1.items():
			if index1 >= index2: sus += index1 * value1
		s2[index2] = sus * 100 / counter
	s1.name = "num"
	s2.name = "cumulative"
	q = s1.to_frame()
	q["cumulative"] = s2
	q.to_csv(output_table)

Blister.EachFile("Clusterization", ["/dev/datasets/FairWind/Pavel/sorted_cut/a*.bam"], "/dev/datasets/FairWind/Pavel/clusterized/", THREADS_NUM = 4)(the_thread)()
