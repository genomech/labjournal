from lib.blister import *

class Test(object):
	
	@staticmethod
	def ClassifiedInput(filemask):
		METHOD_NAME = f"Blister.ClassifiedInput"
		
		input_filenames = Blister.Input([filemask])
		if not input_filenames: return False
		
		length = []
		table = []
		
		for input_filename in input_filenames:
			mods = Blister.FileMods(input_filename)
			table += [mods]
			length += [len(mods)]
		
		length = set(length)
		
		if len(length) > 1:
			print(f"{METHOD_NAME}: Input files cannot be classified.", end='\n')
			return False
		
		cols = list(range(list(length)[0]))
		
		table = pd.DataFrame(table, columns=cols)
		print(table)
		
		sets = []
		reference = 0
		
		for column in cols:
			variables = set(table[column].to_list())
			sets += [len(variables)]
			
		print(sets)

Test.ClassifiedInput("/dev/datasets/ngs_data/biblexome/*.fastq.gz")
