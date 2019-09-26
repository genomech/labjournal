from lib.blister import *

genome = "/dev/datasets/FairWind/_db/rn6.fa"
THREADS_NUM = 12
quality = 30
depth = 20

def the_thread(block, output_dir):
	index, input_filename = block
	
	output_filename = Blister.Output(input_filename, output_dir, f"FilterCalls-D{depth}-Q{quality}", "vcf", rewrite=True)
	if not output_filename: return

	with Blister.Timestamp("FILTER BAM", input_filename, output_filename) as start_time:
				
		command = f"bcftools mpileup -f {genome} {input_filename} | bcftools call -cv -Ou | bcftools filter -i \"DP>{depth} & %QUAL>{quality}\" > {output_filename}"
		#command = f"vcftools --vcf {input_filename} --bed {region} --out {output_filename} --recode"
		sp = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = sp.communicate()
		if err != b'':
			print(f"BamFilter: Shell error: {str(err)}", end='\n')
			return


Blister.EachFile(f"BamFilter", ["/dev/datasets/FairWind/_results/Fatima/sorted/*.bam"], "/dev/datasets/FairWind/_results/Fatima/vcf_D20Q30/", THREADS_NUM)(the_thread)()
