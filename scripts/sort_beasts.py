from lib.blister import *
import vcf



#exit()

# STAGE 2

def the_thread(block, output_dir):
	index, input_filename = block
	return pd.read_csv(input_filename, sep='\t')
results = Blister.EachFile("read_csv", ["/dev/datasets/FairWind/_results/Fatima/comp/*.csv"], "~/", THREADS_NUM = cpu_count())(the_thread)()

table = None

output_filename = Blister.Output("/inter", "/dev/datasets/FairWind/_results/Fatima/comp", "", "csv", rewrite=True)
if not output_filename: exit()

for it in range(len(results)):
	if it == 0:
		table = results[it]
	else:
		table = pd.merge(table, results[it], how='inner', on='CHROM_POS_REF')

table.to_csv(output_filename, sep='\t', index=False)
exit()

# STAGE 1

def the_thread(block, output_dir):
	index, input_filename = block
	lst = []
	output_filename = Blister.Output(input_filename, output_dir, "", "csv", rewrite=True, index=index)
	if not output_filename: return
	with Blister.Read(input_filename, 'r', index) as input_file:
		vcf_reader = vcf.Reader(input_file)
		for record in vcf_reader:
			if record.samples[0].is_het:
				lst += [f"{record.CHROM} {record.POS} {record.REF}"]
		table = pd.DataFrame(lst, columns=['CHROM_POS_REF'])
		table.to_csv(output_filename, sep='\t', index=False)

#mut = ['01', '02', '03', '04', '05', '08']
#mut = [f"/dev/datasets/FairWind/_results/Fatima/vcf_D20Q30/MB_FQ_0{x}_S*_sorted_FilterCalls-D20-Q30.vcf" for x in mut]
#Blister.EachFile("PARSING Mut", mut, "/dev/datasets/FairWind/_results/Fatima/comp/mut", THREADS_NUM = cpu_count())(the_thread)()

wth = ['06', '07', '09', '10']
wth = [f"/dev/datasets/FairWind/_results/Fatima/vcf_D20Q30/MB_FQ_0{x}_S*_sorted_FilterCalls-D20-Q30.vcf" for x in wth]
Blister.EachFile("PARSING WT", wth, "/dev/datasets/FairWind/_results/Fatima/comp/wt_hetero", THREADS_NUM = cpu_count())(the_thread)()
