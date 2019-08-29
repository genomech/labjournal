import pandas as pd
import numpy as np
import time
from multiprocessing import cpu_count, Pool
import functools
import pickle
import matplotlib.pyplot as plt

start_time = time.time()
with open("/dev/datasets/FairWind/_results/bowtie/coverage/merged_table.pickle", 'rb') as f:
    large_table = pickle.load(f)
print(f"Pickle load is done [%f sec]" % (time.time() - start_time), end="\n")

plt.plot(list(large_table[:5000]))
plt.ylabel('Depth average')
plt.xlabel('Distance from restriction site, b')
plt.suptitle('Distance coverage')
plt.savefig("/dev/datasets/FairWind/_results/bowtie/coverage/distance_coverage_5000.svg")

exit()

# --------------------------

def chrom_handler(chrom):
    start_time = time.time()
    
    coverage = pd.read_csv('/dev/datasets/FairWind/_results/bowtie/coverage/grep_coverage/' + chrom + '.txt', header=None, names=['coverage'])
    restrict = pd.read_csv('/dev/datasets/FairWind/_results/bowtie/coverage/restrict/' + chrom + '.txt', header=None, names=['restrict'])
    
    main_table = pd.concat([coverage, restrict], axis=1)
    del coverage
    del restrict
    main_table = main_table.apply(pd.to_numeric, errors='ignore')
    main_table = main_table.groupby(['restrict']).mean()
    
    print(f"{chrom} is done [%.2f sec]" % (time.time() - start_time), end="\n")
    
    return main_table

THREADS_NUM = cpu_count()

chr_list = []
for line in open('/dev/datasets/FairWind/_results/bowtie/coverage/chr_list', 'rt'):
    chr_list += [line[:-1]]

pool = Pool(THREADS_NUM)

results = pool.map(chrom_handler, chr_list)

pool.close()
pool.join()
del pool

start_time = time.time()
large_table = pd.concat(results, axis=1)
del results
print(f"Concat is done [%f sec]" % (time.time() - start_time), end="\n")

# large_table = large_table[:1000000]

start_time = time.time()
large_table = large_table.apply(np.nanmean, axis=1)
print(f"Merging is done [%.2f sec]" % (time.time() - start_time), end="\n")

start_time = time.time()
with open("/dev/datasets/FairWind/_results/bowtie/coverage/merged_table.pickle", 'wb') as f:
    pickle.dump(large_table, f)
print(f"Pickling merged is done [%f sec]" % (time.time() - start_time), end="\n")

start_time = time.time()
large_table.to_csv("/dev/datasets/FairWind/_results/bowtie/coverage/merged_table.csv", sep='\t')
print(f"Pickling merged is done [%f sec]" % (time.time() - start_time), end="\n")

plt.plot(list(large_table[:500]))
plt.ylabel('Depth average')
plt.xlabel('Distance from restriction site, b')
plt.suptitle('Distance coverage')
plt.savefig("/dev/datasets/FairWind/_results/bowtie/coverage/distance_coverage_500.svg")
