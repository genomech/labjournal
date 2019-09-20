from lib.blister import *
import matplotlib.pyplot as plt

def Graph(G_LOG, G_P):
	genes_diff = pd.read_csv("/dev/datasets/FairWind/_results/Fatima/genes_diff.csv", sep='\t')
	genes_diff["sample_1"] = 'wt'
	genes_diff["sample_2"] = 'mut'
	genes_diff["significant"] = genes_diff["significant"].apply(lambda x: True if x == 'yes' else False)
	for col in ['value_1', 'value_2', 'log2_fold_change', 'test_stat', 'p_value', 'q_value']:
		genes_diff[col] = pd.to_numeric(genes_diff[col], downcast='float', errors='raise')
	
	if G_P:
		new_table = genes_diff[genes_diff['q_value'] < 0.05]
		x = new_table["value_1"].to_list()
		y = new_table["value_2"].to_list()
		if G_LOG:
			plt.xscale("log")
			plt.yscale("log")
			axes = plt.gca()
			axes.set_xlim([10,max(x)])
			axes.set_ylim([10,max(y)])
		#plt.rcParams["figure.figsize"] = [12, 9]
		plt.scatter(x, y, s=4)
		hui1 = (", log" if G_LOG else "")
		plt.title(f'Expression in R. norvegicus liver (q < 0.05{hui1})')
		plt.grid(True)
		plt.xlabel('Wild type')
		plt.ylabel('Mutant')
		hui2 = ("_log" if G_LOG else "")
		plt.savefig(f"/dev/datasets/FairWind/_results/Fatima/Expression_all_q005{hui2}.svg")
		plt.clf()
	else:
		sig_true = genes_diff[genes_diff["significant"]]
		sig_false = genes_diff[~genes_diff["significant"]]
		true_list = (sig_true['value_1'].to_list(), sig_true['value_2'].to_list())
		false_list = (sig_false['value_1'].to_list(), sig_false['value_2'].to_list())
		data = (true_list, false_list)
		colors = ("red", "blue")
		groups = ("significant", "non-significant")
		fig = plt.figure()
		ax = fig.add_subplot(1, 1, 1)
		for data, color, group in zip(data, colors, groups):
			x, y = data
			ax.scatter(x, y, c=color, edgecolors='none', s=4, label=group)
			if G_LOG:
				ax.set_yscale('log')
				ax.set_xscale('log')
				axes = plt.gca()
				axes.set_xlim([10,max(x)])
				axes.set_ylim([10,max(y)])
		#plt.rcParams["figure.figsize"] = [12, 9]
		hui3 = (" (log)" if G_LOG else "")
		plt.title(f'Expression in R. norvegicus liver{hui3}')
		plt.grid(True)
		plt.xlabel('Wild type')
		plt.ylabel('Mutant')
		plt.legend(loc=2)
		hui2 = ("_log" if G_LOG else "")
		plt.savefig(f"/dev/datasets/FairWind/_results/Fatima/Expression_all{hui2}.svg")
		plt.clf()
	print("Ready.", end='\n')

Graph(True, True)
Graph(True, False)
Graph(False, True)
Graph(False, False)
