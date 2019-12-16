import pysam
from lib.blister import *

def the_thread(block, output_dir):
	index, input_filename = block
	data = pd.read_csv(input_filename, header=None, names=["1", "2", "3", "4", input_filename], sep="\t")
	return data[input_filename]

results = Blister.EachFile("Concentrator", ["/dev/datasets/FairWind/Pavel/sorted_coverage/animal_*_Coverage.txt"], "/dev/datasets/FairWind/Pavel", THREADS_NUM = cpu_count())(the_thread)()

main = pd.concat(results, axis=1)
main["index"] = pd.Series(list(range(60775505, 60779998)), index=main.index)
main.to_csv("/dev/datasets/FairWind/Pavel/animals_coverage.csv", index=False)
