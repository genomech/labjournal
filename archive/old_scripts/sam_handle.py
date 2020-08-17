from lib.blister import *
import pysam

def the_thread(block, output_dir):
	
	index = block[0]
	item = block[1]
	
	output_filename = Blister.Output(input_filename, output_dir, "UniqueBam", "bam", rewrite=True)
	if not output_filename: return

	with Blister.Read(input_filename, mode='rt', index=item) as input_file, Blister.Timestamp("COPYING", item) as start_time:
		
		samfile = pysam.AlignmentFile(input_file)
		bamfile = pysam.AlignmentFile(output_filename, 'wb', header=samfile.header)
		count = 0
		for read in samfile.fetch():
			bamfile.write(read)
			count += 1
			if count == item * 1000000: break
		
		bamfile.close()

Blister.Logo("Gradient")

input_filenames = Blister.Input(["/dev/datasets/FairWind/_results/bowtie/sam/sample-1-2.sam"])
if not input_filenames: exit()

output_dir = Blister.Dir("/dev/datasets/FairWind/_results/20-120M/bam", create=True)
if not output_dir: exit()

numbers = [20, 40, 60, 80, 100, 120]

with Blister.Threading("GRADIENT") as pool:
	pool.map(functools.partial(the_thread, input_filename=input_filenames[0], output_dir=output_dir), enumerate(numbers))
