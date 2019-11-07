from lib.blister import *
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

fontP = FontProperties()
fontP.set_size('small')
 
table_index = ['Sample',  'Non-coverage, %', 'Middle', 'Median', 'Coverage 75%', 'Coverage 90%', 'Coverage 95%', 'cumulate_list']

def the_thread(block):
	
	index = block[0]
	input_filename = block[1]
	
	with Blister.Timestamp("READ CSV", index) as start_time:
		table = pd.read_csv(input_filename, sep='\t', header=None, names=["chrom", "depth", "a", "b", "percent"])
		table = table.apply(pd.to_numeric, errors='ignore')
		percents_capture = table['percent'].tolist()
		del table
	
	result = pd.Series(index=table_index)
	result['Sample'] = QFileInfo(input_filename).baseName()
	result['Non-coverage, %'] = percents_capture[0] * 100
	result['cumulate_list'] = []
	result['Middle'] = 0
	result['Median'] = 0
	result['Coverage 75%'] = 0
	result['Coverage 90%'] = 0
	result['Coverage 95%'] = 0
	
	with Blister.Timestamp("PROCESSING", index) as start_time:
		for it in range(len(percents_capture)):
			summa = sum(percents_capture[it:])
			result['cumulate_list'] += [summa]
			result['Middle'] += percents_capture[it] * it
			if summa > 0.5: result['Median'] = it
			if summa > 0.75: result['Coverage 75%'] = it
			if summa > 0.9: result['Coverage 90%'] = it
			if summa > 0.95: result['Coverage 95%'] = it
	
	return result

Blister.Logo("Coverage Shorten")

input_filenames = Blister.Input(["/dev/datasets/FairWind/_results/60m/PRIMARY_ANALYSIS_13D/coverage/full/*_ExomeCoverage.txt.all"])
if not input_filenames: exit()

main_table = pd.DataFrame(columns=table_index)

with Blister.Threading("KOVOER") as pool:
	results = pool.map(functools.partial(the_thread), enumerate(input_filenames))
	main_table = main_table.append(results, ignore_index=True)
	del results

with Blister.Timestamp("PLOTTING") as start_time:
	plt.clf()
	for it in main_table.iterrows():
		plt.plot(it[1]["cumulate_list"][:100], label=it[1]["Sample"])
	plt.ylabel('% of A with more depth')
	plt.xlabel('Depth')
	plt.suptitle('Cumulative coverage (60M)')
	plt.legend(prop=fontP)
	plt.savefig("/dev/datasets/FairWind/_results/60m/PRIMARY_ANALYSIS_13D/coverage/full/graph.svg")

main_table.drop(columns=['cumulate_list'], axis=0, inplace=True)
main_table.to_csv("/dev/datasets/FairWind/_results/60m/PRIMARY_ANALYSIS_13D/coverage/full/table.csv")
print(Blister.GitHubTable(main_table))
