from lib.blister import *

genome = "/dev/datasets/FairWind/_db/hg19/hg19.fa"
region = "/dev/datasets/FairWind/_db/MedExome_hg19_capture_targets.bed"
THREADS_NUM = 12

def the_thread(block, output_dir):
	index, input_filename = block
	
	output_filename = Blister.Output(input_filename, output_dir, f"Exome", "vcf", rewrite=True)
	if not output_filename: return

	with Blister.Timestamp("FILTER BAM", input_filename, output_filename) as start_time:
				
		#command = f"bcftools mpileup -f {genome} {input_filename} | bcftools call -cv -Ou | bcftools filter -R {region} -i \"DP>{depth} & %QUAL>{quality}\" > {output_filename}"
		command = f"vcftools --vcf {input_filename} --bed {region} --out {output_filename} --recode"
		sp = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = sp.communicate()
		if err != b'':
			print(f"BamFilter: Shell error: {str(err)}", end='\n')
			return


Blister.EachFile(f"BamFilter", ["/dev/datasets/FairWind/_results/bowtie/vcf_depth10/*.txt"], "/dev/datasets/FairWind/_results/bowtie/vcf_exome_filtered/", THREADS_NUM)(the_thread)()
