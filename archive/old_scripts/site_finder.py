import sys
import math
import subprocess
import re

def processing(chrseq, teststring, filename):
    
    print(f"Start {filename} ... ", end='')
    chrseq = chrseq.upper()
    parts = list(map(len, re.split(teststring, chrseq)));
    if len(parts) <= 1: 
        print("No seqs.", end='\n')
        return False
    
    lst = []
    
    lst += range(int(parts[0]), 0, -1)
    lst += [0] * len(teststring)
    
    for part in parts[1:-1]:
        if part > 1:
            if (part % 2) == 0:
                lst += range(1, int(part / 2) + 1)
                lst += range(int(part / 2), 0, -1)
            else:
                lst += range(1, int(math.floor(part / 2)) + 1)
                lst += range(int(math.ceil(part / 2)), 0, -1)
        elif part == 1:
            lst += [1]
        elif part == 0:
            pass
        
        lst += [0] * len(teststring)
        
    lst += range(1, int(parts[-1]) + 1)
    
    with open(filename, 'wt') as f:
        for item in lst:
            f.write("%d\n" % item)
    
    del lst
    
    print(f"Done.", end='\n')
    

site_r = 'DpnII'
teststring = 'GATC'

genome = 'hg19'
filename = '/dev/datasets/FairWind/_db/hg19/hg19.fa'
outpath = '/dev/datasets/FairWind/_results/bowtie/coverage/restrict/'


chrom = ''
chrseq = ''

for line in open(filename, mode='rt'):

    if line[0] == '>':
        
        processing(chrseq, teststring, outpath + chrom + '.txt')
        
        chrom = line[1:-1]
        chrseq = ''
    else:
        chrseq += line[0:-1]

processing(chrseq, teststring, outpath + chrom + '.txt')
