from lib.blister import *

genome = "/dev/datasets/FairWind/_db/hg19/hg19.fa"
depth = 10
quality = 30
THREADS_NUM = 8

def the_thread(block, output_dir):
	index, input_filename = block
	
	output_filename = Blister.Output(input_filename, output_dir, f"FilterCalls_DP{depth}_QUAL{quality}", "txt", rewrite=True)
	if not output_filename: return

	with Blister.Timestamp("FILTER BAM", input_filename, output_filename) as start_time:
				
		command = f"bcftools mpileup -f {genome} {input_filename} | bcftools call -cv -Ou | bcftools filter -i \"DP>{depth} & %QUAL>{quality}\" > {output_filename}"
		sp = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = sp.communicate()
		if err != b'':
			print(f"BamFilter: Shell error: {str(err)}", end='\n')
			return


Blister.EachFile(f"BamFilter", ["/dev/datasets/FairWind/_results/bowtie/bam/*.bam"], "/dev/datasets/FairWind/_results/bowtie/vcf2/", THREADS_NUM)(the_thread)()
