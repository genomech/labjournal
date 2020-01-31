import pysam
import random

from lib.blister import *

def the_thread(block, output_dir):
	index, input_filename = block
	output_filename = Blister.Output(input_filename, output_dir, "strand-filtered", "bam", rewrite=True, index=index)
	with Blister.Timestamp("Processing", filename_1=input_filename, filename_2=output_filename, index=index):
		samfile = pysam.AlignmentFile(input_filename, "rb")
		filtered = pysam.AlignmentFile(output_filename, "wb", template=samfile)
		current_name = ""
		read_pack = []
		#total = 0
		for read in samfile:
			if read.is_unmapped:
				filtered.write(read)
				continue
			if read.query_name != current_name:
				for val in set([str(x.reference_start) + ":" + str(x.reference_length) for x in read_pack]):
					l = []
					for x in [(str(x.reference_start) + ":" + str(x.reference_length), x) for x in read_pack]:
						if x[0] == val: l += [x[1]]
					filtered.write(random.choice(l))
				current_name = read.query_name
				read_pack = []
			
			read_pack += [read]
			#total += 1
			#if total == 100: break
		samfile.close()
		filtered.close()
	
Blister.EachFile("Кусь за Русь!", ["/dev/datasets/FairWind/_results/case/128_S5_dupless.bam"], "/dev/datasets/FairWind/_results/case/", THREADS_NUM = cpu_count())(the_thread)()
