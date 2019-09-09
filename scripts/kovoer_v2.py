from blister import *
import pandas as pd
from tabulate import tabulate
import matplotlib.pyplot as plt
 
table_index = ['Sample',  'Non-coverage, %', 'Middle', 'Median', 'Coverage 75%', 'Coverage 90%', 'Coverage 95%', 'cumulate_list']

def the_thread(block):
	
	index = block[0]
	input_filename = block[1]
	
	with blister_timestamp("READ CSV", index) as start_time:
		table = pd.read_csv(input_filename, sep='\t')
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
	
	with blister_timestamp("PROCESSING", index) as start_time:
		for it in range(len(percents_capture)):
			summa = sum(percents_capture[it:])
			result['cumulate_list'] += [summa]
			result['Middle'] += percents_capture[it] * it
			if summa > 0.5: result['Median'] = it
			if summa > 0.75: result['Coverage 75%'] = it
			if summa > 0.9: result['Coverage 90%'] = it
			if summa > 0.95: result['Coverage 95%'] = it
	
	return result

blister_logo("Coverage Shorten 60M")

input_filenames = blister_input(["/dev/datasets/FairWind/_results/60m/shorten_tables/*.csv"])
if not input_filenames: exit()

main_table = pd.DataFrame(columns=table_index)

with blister_threading("KOVOER") as pool:
	results = pool.map(functools.partial(the_thread), enumerate(input_filenames))
	main_table = main_table.append(results, ignore_index=True)
	del results

with blister_timestamp("PLOTTING") as start_time:
	plt.clf()
	for it in main_table.iterrows():
		plt.plot(it[1]["cumulate_list"][:100], label=it[1]["Sample"])
	plt.ylabel('% of A with more depth')
	plt.xlabel('Depth')
	plt.suptitle('Cumulative coverage (all)')
	plt.legend()
	plt.savefig("/dev/datasets/FairWind/_results/60m/allplot_60m_fq.svg")

main_table.drop(columns=['cumulate_list'], axis=0, inplace=True)
print(tabulate(main_table, headers=main_table.columns, tablefmt="github", showindex=False))
