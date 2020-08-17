import pysam
from lib.blister import *

def the_thread(block, output_dir):
	index, input_filename = block
	output_filename = Blister.Output(input_filename, output_dir, "", "bam", rewrite=True, index=index)
	samfile = pysam.AlignmentFile(input_filename, "r")
	pairedreads = pysam.AlignmentFile(output_filename, "wb", template=samfile)
	for read in samfile:
		if read.query_name.find("__") != -1: read.query_name = read.query_name[:-10]
		#print(read)
		pairedreads.write(read)

results = Blister.EachFile("TrashCutter", ["/dev/datasets/FairWind/Pavel/sorted/_animal_0.sam"], "/dev/datasets/FairWind/Pavel/sorted_cut", THREADS_NUM = cpu_count())(the_thread)()
