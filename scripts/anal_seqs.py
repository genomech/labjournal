import gzip
import sys
import pandas as pd
import numpy as np
import string
import re

def Out(total_):
    print("Total: %d" % (total_), end='\r')


filename = '/dev/datasets/FairWind/_results/35/sample-1-5_R1_Illuminaless.fastq'
output_file = '/dev/datasets/FairWind/_results/35/sample-1-5_R1_Illuminaless_STATISTICS.txt'
seqs = pd.DataFrame(np.array([['-bridge-', 'GCTGAGG'], ['-egdirb-', 'CCTCAGC'], ['-gatc-', 'GATC']]), columns=['name', 'seq'])
genome = '-genome-'

print(f"\nHi there.\nWe're mapping reads with:\n")
print(seqs)

input0 = open(filename, 'r')

counter = 0
total = 0
tyk = 0
main_list = []
main_table = pd.DataFrame(columns=['count', 'mask'])

pd.set_option('max_colwidth',1000)

for it in range(4000000):
    
    line = input0.readline()

    Out(total)

    counter += 1
    if counter == 5:
        counter = 1

    if (counter == 2):
        line = line[0:-1]
        for index, row in seqs.iterrows():
            line = line.replace(row['seq'], row['name'])
        for it in range(10):
            its = str(it + 1)
            line = re.sub("(^|-)[ATGCN]{" + its + "}($|-)", "--[" + its + "]--", line)

        line = re.sub(r"[ATGCN]{11,}", genome, line)
        # optimization
        line = line.replace('ome--gatc--gen', '')
        line = re.sub(r"^-gatc--genome-", genome, line)
        line = re.sub(r"-genome--gatc-$", genome, line)

        main_list.append(line)
        total += 1

shorted_list = list(set(main_list))
for note in shorted_list:
    main_table = main_table.append(pd.Series([main_list.count(note), note], index=['count', 'mask']), ignore_index=True)

output0 = open(output_file, 'w')
main_table.sort_values(by=['count'], ascending=False).to_string(output0)
print('\n')

input0.close()
output0.close()
