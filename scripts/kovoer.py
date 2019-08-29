import pickle
import pandas as pd
import matplotlib.pyplot as plt

with open("../_pickles/csv/MedExome_coverage_withdup_gist_[all].pd.pickle", 'rb') as f1:
    hist_capture = pickle.load(f1)

hist_capture = hist_capture.apply(pd.to_numeric, errors='ignore')
print (hist_capture)

#with open("../_pickles/csv/hist_hg19_[all].pd.pickle", 'rb') as f2:
#    hist_hg19 = pickle.load(f2)

#hist_hg19 = hist_hg19.apply(pd.to_numeric, errors='ignore')

print("Pickle loaded.", end="\n")

percents_capture = hist_capture['% of A at depth'].tolist()
#percents_hg19 = hist_hg19['% of A at depth'].tolist()


plt.plot(percents_capture[:100])
plt.ylabel('% from A at depth')
plt.xlabel('Depth')
plt.suptitle('Coverage of capture')
plt.savefig("../_graph/coverage_100_capture_wd.svg")

plt.clf()

#plt.plot(percents_hg19[:100])
#plt.ylabel('% from A at depth')
#plt.xlabel('Depth')
#plt.suptitle('Coverage of hg19')
#plt.savefig("../_graph/coverage_100_hg19.svg")

print("Coverage is ready.", end="\n")

cover_50_capture = 0
cover_75_capture = 0
cover_90_capture = 0
cover_95_capture = 0

#cover_50_hg19 = 0
#cover_75_hg19 = 0
#cover_90_hg19 = 0
#cover_95_hg19 = 0

cumulate_capture = []
#cumulate_hg19 = []
middle_death_capture = 0
#middle_death_hg19 = 0

for it in range(len(percents_capture)):
    summa = sum(percents_capture[it:])
    cumulate_capture += [summa]
    middle_death_capture += percents_capture[it] * it
    if summa > 0.5: cover_50_capture = it
    if summa > 0.75: cover_75_capture = it
    if summa > 0.9: cover_90_capture = it
    if summa > 0.95: cover_95_capture = it

#for it in range(len(percents_hg19)):
#    summa = sum(percents_hg19[it:])
#    cumulate_hg19 += [summa]
#    middle_death_hg19 += percents_hg19[it] * it
#    if summa > 0.5: cover_50_hg19 = it
#    if summa > 0.75: cover_75_hg19 = it
#    if summa > 0.9: cover_90_hg19 = it
#    if summa > 0.95: cover_95_hg19 = it

plt.clf()

plt.plot(cumulate_capture[:100])
#plt.plot(cumulate_hg19[:100])
plt.ylabel('% of A with more depth')
plt.xlabel('Depth')
plt.suptitle('Cumulative coverage')
plt.savefig("../_graph/cumulate_100_wd.svg")

print("Cumulate is ready.", end="\n")

print("\nCAPTURE DATA", end="\n")

print(f"Non-coverage = {percents_capture[0] * 100}%", end="\n")
print(f"Median = depth {cover_50_capture}", end="\n")
print(f"Middle = depth {middle_death_capture}", end="\n")
print(f"Cover 75% = depth {cover_75_capture}", end="\n")
print(f"Cover 90% = depth {cover_90_capture}", end="\n")
print(f"Cover 95% = depth {cover_95_capture}", end="\n")

#print("\nHG19 DATA", end="\n")
#
#print(f"Non-coverage = {percents_hg19[0] * 100}%", end="\n")
#print(f"Median = depth {cover_50_hg19}", end="\n")
#print(f"Middle = depth {middle_death_hg19}", end="\n")
#print(f"Cover 75% = depth {cover_75_hg19}", end="\n")
#print(f"Cover 90% = depth {cover_90_hg19}", end="\n")
#print(f"Cover 95% = depth {cover_95_hg19}", end="\n")
