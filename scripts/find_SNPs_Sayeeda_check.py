from lib.blister import *
import vcf

checker_path = "/dev/datasets/FairWind/_results/Fatima/checkers/1.csv"
mut = ['01', '02', '03', '04', '05', '08']
wth = ['06', '07', '09', '10']

def merger(table_list):

	table = None
	for it in range(len(table_list)):
		if it == 0:
			table = table_list[it]
		else:
			table = pd.merge(table, table_list[it], how='inner', on='CHROM_POS_REF')
	return table

def the_thread(block, output_dir):
	index, input_filename = block
	lst = []
	output_filename = Blister.Output(input_filename, output_dir, "", "csv", rewrite=True, index=index)
	if not output_filename: return
	with Blister.Read(input_filename, 'r', index) as input_file:
		vcf_reader = vcf.Reader(input_file)
		for record in vcf_reader:
			if output_dir[-2:] == '11':
				if not record.samples[0].is_het:
					lst += [f"{record.CHROM} {record.POS} {record.REF}"]
			elif output_dir[-2:] == '10':
				if record.samples[0].is_het:
					lst += [f"{record.CHROM} {record.POS} {record.REF}"]
		table = pd.DataFrame(lst, columns=['CHROM_POS_REF'])
		return table

mut = [f"/dev/datasets/FairWind/_results/Fatima/vcf_D20Q30/MB_FQ_0{x}_S*_sorted_FilterCalls-D20-Q30.vcf" for x in mut]
mut_results = Blister.EachFile("PARSING Mut", mut, "11", THREADS_NUM = cpu_count())(the_thread)()

wth = [f"/dev/datasets/FairWind/_results/Fatima/vcf_D20Q30/MB_FQ_0{x}_S*_sorted_FilterCalls-D20-Q30.vcf" for x in wth]
wth_results = Blister.EachFile("PARSING WT", wth, "10", THREADS_NUM = cpu_count())(the_thread)()

mut_table = merger(mut_results)
wth_table = merger(wth_results)

snps = pd.merge(mut_table, wth_table, how='inner', on='CHROM_POS_REF')
snps.to_csv(checker_path, sep='\t', index=False)
