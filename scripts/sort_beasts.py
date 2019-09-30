from lib.blister import *
import vcf

Blister.Logo("Find Beasts")

all_wt = pd.read_csv('/dev/datasets/FairWind/_results/Fatima/comp/all_wt.csv', sep='\t')
homo_mut = pd.read_csv('/dev/datasets/FairWind/_results/Fatima/comp/homo_mut.csv', sep='\t')

great_table = pd.DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT'])

count = 0
with Blister.Timestamp("DO") as start_time:
	homo_mut = homo_mut[~all_wt.isin(homo_mut)]['POS'].dropna().to_list()

print(homo_mut)

exit()

i = ['/dev/datasets/FairWind/_results/Fatima/comp/1.csv', '/dev/datasets/FairWind/_results/Fatima/comp/2.csv']

input_filenames = Blister.Input(i)
if not input_filenames: exit()

output_dir = Blister.Dir('/dev/datasets/FairWind/_results/Fatima/comp/')
if not output_dir: exit()

output_filename = Blister.Output("/3", output_dir, "", "csv", rewrite=True)
if not output_filename: exit()

lst = []

for input_filename in input_filenames:
	lst += [pd.read_csv(input_filename, sep='\t')]

great_table = None

for it in range(len(lst)):
	if it == 0: great_table = lst[it]
	else: great_table = pd.merge(great_table, lst[it], how='inner', on=['CHROM', 'POS', 'REF', 'ALT'])

great_table.sort_values(by=['CHROM', 'POS'], inplace=True)
great_table.to_csv(output_filename, sep='\t', index=False)
