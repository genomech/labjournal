from sklearn.cluster import AgglomerativeClustering as AC
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist
from sklearn import metrics
from scipy.cluster.hierarchy import fcluster
from scipy import stats
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import Levenshtein
import json
import pandas

with open("/dev/datasets/FairWind/.cloud/FORGE/data_file.json", "r") as write_file: Connections = json.load(write_file)
check = list(Connections)[:22]
table = pandas.DataFrame(index=check, columns=check)

for i in table.columns:
	for j in table.index:
		table.loc[j,i] = Levenshtein.distance(j,i)
table.to_csv("/dev/datasets/FairWind/.cloud/FORGE/Levenstein.csv")
TopicData = table.values.tolist()

#table["Ones"] = table.apply(lambda x: x.to_list().count(1), axis=1)
table["Ones_list"] = table.apply(lambda x: [index for index, item in enumerate(x.to_list()) if item <= 2], axis=1)
print(table["Ones_list"].to_list())

TopicData = pdist(TopicData)
#print(TopicData)
Labels = [str(index) + "(" + str(table.loc[item, "Ones"]) + "): " + item for index, item in enumerate(table.index.to_list())]

LinkageMatrix = hierarchy.linkage(TopicData, method="single", metric="euclidean")

plt.figure(figsize=(12,12))
plt.title("Seqs clustered", fontsize=14)
plt.xlabel("Distance (ward/Euclidean)", fontsize=12)
matplotlib.rcParams['lines.linewidth'] = 1.2
hierarchy.dendrogram(LinkageMatrix, truncate_mode="level", color_threshold=8.0, show_leaf_counts=True, no_labels=False, labels=Labels, orientation="left")
plt.savefig("/dev/datasets/FairWind/.cloud/FORGE/cluster_analysis.svg", dpi=300, figsize=(12,18), bbox_inches="tight")
