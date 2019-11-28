from Bio import SeqIO
from lib.blister import *
import pysam

Blister.Logo("Genome Surfer s2 v0.01")

def the_thread(block, output_dir):
	index, input_filename = block
	main_table = pd.read_csv(input_filename, sep='\t')
	
	gatcF = main_table[(main_table["Strand"] == True) & (main_table["Context"] == "GATC")].shape[0] * 100 / main_table[(main_table["Strand"] == True)].shape[0]
	gatcR = main_table[(main_table["Strand"] == False) & (main_table["Context"] == "GATC")].shape[0] * 100 / main_table[(main_table["Strand"] == False)].shape[0]
	
	posF = pd.DataFrame(columns=['A', 'T', 'G', 'C', 'N'], index=[1, 2, 3, 4])
	posR = pd.DataFrame(columns=['A', 'T', 'G', 'C', 'N'], index=[1, 2, 3, 4])
	posF.fillna(0, inplace=True)
	posR.fillna(0, inplace=True)
	
	for i in main_table.iterrows():
		for k in [1, 2, 3, 4]:
			if i[1]["Context"] != i[1]["Context"]:
				posF['N'][k] += 1
				continue
			if i[1]["Strand"]: posF[i[1]["Context"][k - 1]][k] += 1
			else: posR[i[1]["Context"][k - 1]][k] += 1
	
	posF = posF.apply(lambda x: x * 100 / main_table[(main_table["Strand"] == True)].shape[0])
	posR = posR.apply(lambda x: x * 100 / main_table[(main_table["Strand"] == False)].shape[0])
	print(f"{input_filename}. GATC: F = {gatcF}, R = {gatcR}\nF:\n{posF.to_string()}\nR:\n{posR.to_string()}")
	
Blister.EachFile("Processing", ["/dev/datasets/FairWind/_results/November/November_BGG_bam_context/191107_X603_FCH5KNCCCX2_L5_*_BGG-100k_sorted_context.csv"], "/dev/datasets/FairWind/_results/November/November_BGG_bam_context")(the_thread)()

exit()

# ---------------------------

Blister.Logo("Genome Surfer v0.01")

genome = {}

with Blister.Timestamp("GENOME 2 DICT") as start_time:
	for record in SeqIO.parse("/dev/datasets/FairWind/_db/hg19/hg19.fa", "fasta"): genome[record.id] = record.seq.__str__()

def the_thread(block, output_dir):
		index, input_filename = block
		output_filename = Blister.Output(input_filename, output_dir, "context", "csv", rewrite=True, index=index)
		main_table = pd.DataFrame(columns=["ID", "Strand", "Context"])
		samfile = pysam.AlignmentFile(input_filename, "rb")
		for read in samfile:
			if read.reference_name is None: continue
			context = genome[read.reference_name][read.reference_end : read.reference_end + 4] if read.is_reverse else genome[read.reference_name][read.reference_start - 4 : read.reference_start]
			ser = pd.Series([read.query_name, (not read.is_reverse), context.upper()], index=["ID", "Strand", "Context"])
			main_table = main_table.append(ser, ignore_index=True)
		samfile.close()
		main_table.to_csv(output_filename, sep='\t', index=False)

Blister.EachFile("Processing", ["/dev/datasets/FairWind/_results/November/November_BGG_bam/*.bam"], "/dev/datasets/FairWind/_results/November/November_BGG_bam_context")(the_thread)()
