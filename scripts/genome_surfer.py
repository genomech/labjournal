from Bio import SeqIO
from lib.blister import *
import pysam

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
