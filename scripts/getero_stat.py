from lib.blister import *
import vcf
import numpy as np
import matplotlib.pyplot as plt

filenames = ['/dev/datasets/FairWind/_results/bowtie/duplicates_experiment/sample-1-5_*_exome.vcf.recode.vcf']
dir_path = '/dev/datasets/FairWind/_results/bowtie/duplicates_experiment/'

def lister(table, less):
	if less:
		table = table[table['DP'] < 50]
	else:
		table = table[table['DP'] >= 50]
	table['lol'] = 1
	table.drop(columns=['DP', 'REF-DP', 'ALT-DP'], inplace=True)
	table2 = table.groupby(by=['ALT-F']).apply(np.sum)
	table2.drop(columns=['ALT-F'], inplace=True)
	table2 = table2.to_dict()['lol']
	for i in range(100):
		if not (i in table2):
			table2[i] = 0
	return pd.Series(table2).sort_index().apply(int).to_list()

def the_thread(block, output_dir):
	index, input_filename = block
	
	output_filename_1 = Blister.Output(input_filename, output_dir, "10-50", "svg", rewrite=True, index=index)
	if not output_filename_1: return

	output_filename_2 = Blister.Output(input_filename, output_dir, "50-inf", "svg", rewrite=True, index=index)
	if not output_filename_2: return
	
	table = pd.DataFrame(columns=['DP', 'REF-DP', 'ALT-DP', 'ALT-F'])
	
	with Blister.Timestamp("Parsing VCF", input_filename, index=index) as start_time, Blister.Read(input_filename, mode='r', index=index) as input_file:
		vcf_reader = vcf.Reader(input_file)
		total = 0
		for record in vcf_reader:
			total += 1
			#if total > 1000: break
			if record.samples[0].is_het:
				dp4 = record.INFO['DP4']
				ref = dp4[0] + dp4[1]
				alt = dp4[2] + dp4[3]
				table = table.append(pd.Series([record.INFO['DP'], ref, alt, int(alt / record.INFO['DP'] * 100)], index=['DP', 'REF-DP', 'ALT-DP', 'ALT-F']), ignore_index=True)
	
		table_less = lister(table, True)
		table_greater = lister(table, False)
		
		return [output_filename_1, output_filename_2, table_less, table_greater]
	
results = Blister.EachFile("Hetero Statistics", filenames, dir_path)(the_thread)()

with Blister.Timestamp("Plotting") as start_time:
		name = results[0][0].split('/')[-1].split('.')[0]
		plt.clf()
		plt.plot(results[0][2], label=results[0][0].split('/')[-1].split('.')[0])
		plt.plot(results[1][2], label=results[1][0].split('/')[-1].split('.')[0])
		plt.ylabel('Number of calls')
		plt.xlabel('% alt allele')
		plt.suptitle(f'Heterozygote stat (10..50)')
		plt.legend()
		plt.savefig(results[0][0])
		
		name = results[0][1].split('/')[-1].split('.')[0]
		plt.clf()
		plt.plot(results[0][3], label=results[0][1].split('/')[-1].split('.')[0])
		plt.plot(results[1][3], label=results[1][1].split('/')[-1].split('.')[0])
		plt.ylabel('Number of calls')
		plt.xlabel('% alt allele')
		plt.suptitle(f'Heterozygote stat (50+)')
		plt.legend()
		plt.savefig(results[0][1])
