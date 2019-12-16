import vcf
from lib.blister import *

def the_thread(block, output_dir):
	index, input_filename = block
	output_filename = Blister.Output(input_filename, output_dir, "QUAL20", "vcf", rewrite=True, index=index)
	with Blister.Timestamp("Filtering", filename_1=input_filename, filename_2=output_filename, index=index):
		vcf_reader = vcf.Reader(open(input_filename, 'r'))
		vcf_writer = vcf.Writer(open(output_filename, 'w'), template=vcf_reader)
		for record in vcf_reader:
			if record.QUAL > 20: vcf_writer.write_record(record)

Blister.EachFile("Filter QUAL", ["/dev/datasets/FairWind/_results/cut/filtercalls/*.vcf"], "/dev/datasets/FairWind/_results/cut/filtercalls_qual20", THREADS_NUM = cpu_count())(the_thread)()
