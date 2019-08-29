import pickle
import pandas as pd

with open("../_pickles/csv/Ensembl_transcripts_[all].pd.pickle", 'rb') as f:
    transcripts = pickle.load(f)
transcripts['Transcript end (bp)'] = transcripts['Transcript end (bp)'].apply(pd.to_numeric, errors='ignore')
transcripts['Transcript start (bp)'] = transcripts['Transcript start (bp)'].apply(pd.to_numeric, errors='ignore')

transcripts['length'] = transcripts['Transcript end (bp)'] - transcripts['Transcript start (bp)']

new_transcripts = pd.DataFrame(columns=transcripts.columns)

dick = transcripts['Gene stable ID'].unique()

did = 0
total = len(dick)

for cock in dick:
    table = transcripts[transcripts['Gene stable ID'] == cock]
    new_transcripts = new_transcripts.append(table.loc[table['length'].idxmax()], ignore_index=True)
    did += 1
    #if did == 20: break
    print("%.2f%%" % (did * 100 / total), end='\r')

print("", end='\n')

#with open("../_pickles/csv/Ensembl_transcripts_[good].pd.pickle", 'wb') as f:
#    pickle.dump(new_transcripts, f)

new_transcripts.to_csv("./Ensembl_maxlength.csv", sep="\t", index=False)

print("Done", end='\n')
print(new_transcripts)
