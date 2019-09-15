from blister import *

filenames = ["/dev/datasets/FairWind/_results/bowtie/bam/*.bam"]
dir_path = "/dev/datasets/FairWind/_results/bowtie/vcf/"
genome = "/dev/datasets/FairWind/_db/hg19/hg19.fa"
depth = 10
quality = 30
THREADS_NUM = cpu_count()

Blister.Logo(f"Monster BamFilter")
		
input_filenames = Blister.Input(filenames)
if not input_filenames: exit() 
		
output_dir = Blister.Dir(dir_path, create=True)
if not output_dir: exit()
		
for input_filename in input_filenames:
			
	output_filename = Blister.Output(input_filename, output_dir, f"FilterCalls_DP{depth}_QUAL{quality}", "txt", rewrite=True)
	if not output_filename: continue

	with Blister.Timestamp("FILTER BAM", input_filename, output_filename) as start_time:
				
		command = f"bcftools mpileup --threads {THREADS_NUM} -f {genome} {input_filename} | bcftools call --threads {THREADS_NUM} -cv -Ou | bcftools filter -i \"DP>{depth} & %QUAL>{quality}\" > {output_filename}"
		sp = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = sp.communicate()
		if err != b'':
			print(f"{METHOD_NAME}: Shell error: {str(err)}", end='\n')
			continue
