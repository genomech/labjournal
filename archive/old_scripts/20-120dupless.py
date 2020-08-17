import pysam

from lib.blister import *

lens = [20000000, 40000000, 60000000, 80000000, 100000000, 120000000]
out_dir = Blister.Dir("/dev/datasets/FairWind/_results/cut/uncut_picard/20-120strandless")
input_filename = "/dev/datasets/FairWind/_results/cut/uncut_picard/strandless/sample-1-5_strand-filtered.bam"

def the_thread(block):
		index, item = block
		output_filename = Blister.Output(input_filename, out_dir, str(item), "bam", rewrite=True, index=index)
		total = 0
		with Blister.Timestamp("Processing", filename_1=output_filename, index=index):
			samfile = pysam.AlignmentFile(input_filename, "rb")
			filtered = pysam.AlignmentFile(output_filename, "wb", template=samfile)
			for read in samfile:
				filtered.write(read)
				total += 1
				if total == item: break
			samfile.close()
			filtered.close()

Blister.Logo("20-120 Strandless")

with Blister.Threading("20-120 Strandless", THREADS_NUM = cpu_count()) as pool:
	pool.map(functools.partial(the_thread), enumerate(lens))
