from pybedtools import BedTool

names = BedTool('/dev/datasets/FairWind/_db/rn6.bed')
nicks = BedTool('/dev/datasets/FairWind/_results/Fatima/cuffmerge/merged.gtf')

a_and_b = nicks.intersect(names, names=True)
print(a_and_b)
