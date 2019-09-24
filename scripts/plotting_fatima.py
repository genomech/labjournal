from lib.blister import *
import matplotlib.pyplot as plt
import numpy as np

labels = ['006_S18[wt]', '007_S19[wt]', '009_S21[wt]', '010_S22[wt]', '001_S13[mut]', '002_S14[mut]', '003_S15[mut]', '004_S16[mut]', '005_S17[mut]', '008_S20[mut]']
# values = [3.80141, 18.945, 9.62214, 7.51061, 652.36, 573.313, 286.91, 215.396, 466.817, 3927.61] # hdc
#values = [0, 0, 0, 0, 1907.49, 0, 0, 0, 0, 0] # magi2
#values = [0.151139, 0.185005, 0.0620009, 0.0906103, 16.4181, 24.9628, 23.1145, 29.6679, 53.4047, 4.7837] # fgf23
#values = [0.379374, 0.544733, 0.899505, 0.590493, 12.7651, 15.2278, 17.4665, 18.6676, 20.7459, 132.254] # tnni1
values = [0.871171, 0.702152, 1.14676, 0.776064, 13.848, 24.3044, 6.70644, 11.1848, 57.7782, 94.6212] # atf3
x = np.arange(len(labels))  # the label locations

plt.rcdefaults()
fig, ax = plt.subplots()

ax.barh(x, values, align='center')
ax.set_yticks(x)
ax.set_yticklabels(labels)
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel('FPKM')
ax.set_title('Expression in Atf3 gene')
plt.savefig("/dev/datasets/FairWind/.cloud/core/labjournal/scripts_results/sayeeda_atf3.svg")

exit()

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
