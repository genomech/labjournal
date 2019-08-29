import pandas as pd
from multiprocessing import cpu_count, Pool
import functools
import time
import pickle

in_path = '/dev/datasets/FairWind/_results/dinara/new_pickles/'
filenames = ['le1_S7', 'le2_S8'] # .pickle
out_path = '/dev/datasets/FairWind/_results/dinara/csv_merged_light/'

print(f"\n=== MARGO 0.1 ===\n", end="\n")
start_time = time.time()

with open(f"{in_path}{filenames[0]}_light.pickle", 'rb') as f1:
    table1 = pickle.load(f1)
    
print(f"Pickle {filenames[0]} is loaded [%f sec]" % (time.time() - start_time), end="\n")
start_time = time.time()

with open(f"{in_path}{filenames[1]}_light.pickle", 'rb') as f2:
    table2 = pickle.load(f2)
    
print(f"Pickle {filenames[1]} is loaded [%f sec]" % (time.time() - start_time), end="\n")
start_time = time.time()

table_merged = pd.merge(table1, table2, how='inner')

del table1
del table2

print(f"Merging tables is done [%f sec]" % (time.time() - start_time), end="\n")
start_time = time.time()

table_merged.to_csv(out_path + filenames[0] + "-" + filenames[1] + "_merged_light.csv", sep='\t', index=False, mode='w')
    
print(f"Writing to CSV file is done [%f sec]" % (time.time() - start_time), end="\n")
