from lib.Liebe import *
import pandas as pd
import re
import pickle

def parse_header(info):
	re_list = re.split("\|", info)
	re_list[0] = re.findall('\S+$', re_list[0])[0]
	re_list[-1] = re.findall('^[^"]+', re_list[-1])[0]
	return re_list

def parse_line(line, vcf_header, csq_header):
	line_list = re.split("\t", line)
	info = line_list[7]
	sample = line_list[9]
	del line_list[slice(-3,-1)]
	del line_list[-1]
	line_list += re.split(":", sample)
	line_list[-1] = (line_list[-1])[0:-1]
	gt_list = re.split("/", line_list[-2])
	if len(gt_list) == 2:
		if gt_list[0] == gt_list[1]:
			line_list[-2] += " [homo]"
		else:
			line_list[-2] += " [hetero]"
    
	chrom, pos = line_list[0], line_list[1]

	csq_block = pd.DataFrame(columns=csq_header)
	csq_bigline = re.split("CSQ=", info)[1]
	csq_lines = re.split(",", csq_bigline)
	for csq_line in csq_lines:
		csq_line_list = re.split("\|", csq_line)
		if len(csq_line_list) == len(csq_header):
			csq_block = csq_block.append(pd.Series(csq_line_list, index=csq_header), ignore_index=True)
		else:
			csq_block = csq_block.append(pd.Series(['DAMAGED'] * len(csq_header), index=csq_header), ignore_index=True)
			print(' - DAMAGED string found:\nFull: %s\nPart: %s', end="\n")
        
	vcf_block = pd.DataFrame(columns=vcf_header)
	for it in range(len(csq_block.index)):
		vcf_block = vcf_block.append(pd.Series(line_list, index=vcf_header), ignore_index=True)
    
	table = pd.concat([vcf_block, csq_block], axis=1)
	return table

# __main__

Blister.Logo("Blistered VAPE")

input_filenames = Blister.Input(["/dev/datasets/FairWind/_results/bowtie/sample-1-9_NF1-2_vep.vcf"])
if not input_filenames: exit()

hgmd_filename = Blister.Input(["/dev/datasets/FairWind/_db/hgmd.db"])[0]

pickles_output_dir = Blister.Dir("/dev/datasets/FairWind/_results/bowtie/", create=True)
if not pickles_output_dir: exit()
tables_output_dir = Blister.Dir("/dev/datasets/FairWind/_results/bowtie/", create=True)
if not tables_output_dir: exit()

with Blister.Timestamp("HGMD table") as start_time:

	hgmd_data = pd.read_csv(hgmd_filename, sep=',')
	hgmd_data.rename(columns={"Variant name":"HGMD", "Chromosome/scaffold name":"CHROM", "Chromosome/scaffold position start (bp)":"POS"}, inplace=True)
	hgmd_data["CHROM"] = hgmd_data["CHROM"].apply(lambda x: "chr"+str(x))
	hgmd_data.drop(columns={'Chromosome/scaffold position end (bp)', 'Variant source'}, axis=1, inplace=True) 
	hgmd_data["POS"] = hgmd_data["POS"].apply(lambda x: str(x))

for input_filename in input_filenames:
	
	pickle_filename = Blister.Output(input_filename, pickles_output_dir, "", "pickle", rewrite=True)
	if not pickle_filename: continue
	table_filename = Blister.Output(input_filename, tables_output_dir, "", "csv", rewrite=True)
	if not table_filename: continue

	with Blister.Read(input_filename, 'rt') as input_file:
		
		with Blister.Timestamp("Header") as start_time:
			csq_header = list()
			vcf_header = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'GT', 'PL']
			while 1:
				line = input_file.readline()
				if (re.match("^#[^#].*$", line) != None):
					break
				if (re.match("^##INFO=<ID=CSQ.*$", line) != None):
					csq_header = parse_header(line)
			main_table = pd.DataFrame(columns=(vcf_header + csq_header))
    
		with Blister.Threading("Parsing VCF") as pool:
			results = pool.map(functools.partial(parse_line, vcf_header=vcf_header, csq_header=csq_header), input_file)
			main_table = main_table.append(results, ignore_index=True)
			del results

		with Blister.Timestamp("Merge with HGMD") as start_time:
			main_table = pd.merge(hgmd_data, main_table, how='right', on=["CHROM","POS"])
			main_table.sort_values(by=["CHROM","POS"], inplace=True)
    
		with Blister.Timestamp("Pickling", filename_1=pickle_filename) as start_time, open(pickle_filename, 'wb') as f:
			pickle.dump(main_table, f)

		with Blister.Timestamp("Save table", filename_1=input_filename, filename_2=table_filename) as start_time:
			main_table.to_csv(table_filename, sep='\t', index=False, mode='w')
			del csq_header
			del vcf_header
			del main_table
