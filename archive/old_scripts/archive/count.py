import sys

def Out(found_, total_):
    print("Found: %d | Total: %d" % (found_, total_), end='\r')

filename = '/dev/datasets/FairWind/_results/35/sample-1-3_R1_Illuminaless_TO.fastq'

input0 = open(filename, 'r')

new_found_land = ""
counter = 0
total = 0
found = 0

for t in range(4000000):

    line = input0.readline()
    
    Out(found, total)

    counter += 1
    if counter == 5:
        counter = 1
    if counter < 1:
        continue

    tyk = line.find("CTCAGCGCTGAG")

    if (counter == 2) and (tyk == -1):
        new_found_land = ""
        counter = -2
        total += 1
        continue

    if counter == 4:
        new_found_land = ""
        found += 1
        total += 1

print('\n')

input0.close()
output0.close()
