import gzip
import sys
import numpy as np
import string

def Out(total_):
    print("Total: %d" % (total_), end='\r')

filename = '/dev/datasets/ngs_data/Chinese/bat-Hi-C_mES_alu1_1_10M.fastq'
seq = "CGGTGGC"
length = len(seq)
poses = list([-12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, length + 1, length + 2, length + 3, length + 4, length + 5, length + 6, length + 7, length + 8, length + 9, length + 10, length + 11, length + 12])
nuc = list('ATGCN0')
a = np.zeros(shape=(150,24,6))

print(f"\nPalindrome deep analysis.\nMaking a dump...\n")

input0 = open(filename, 'r')

counter = 0
total = 0
tyk = 0

for line in input0:

    Out(total)

    counter += 1
    if counter == 5:
        counter = 1

    if (counter == 2):
        line_d = line
        tyk = line_d.find(seq)
        if not (tyk == -1):
            for i in range(len(poses)):
                pos_abs = tyk + poses[i]
                if (pos_abs < 0) or (pos_abs > 149):
                    a[tyk, i, nuc.index('0')] += 1
                else:
                    a[tyk, i, nuc.index(line_d[pos_abs])] += 1

        total += 1

print(np.sum(a, axis = (1, 2)))
np.save('/dev/datasets/FairWind/_results/chinese.dump', a)

print('\n')

input0.close()
