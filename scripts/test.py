from lib.blister import *

class Huita(object):
	
	@staticmethod
	def ClassifyInput(filenames, sep='_'):
		if not filenames: return False
		METHOD_NAME = f"Blister.ClassifyInput"
		lls = []
		for filename in filenames:
			fileinfo = QFileInfo(filename)
			dir_path, basename, suffix = [fileinfo.absolutePath(), fileinfo.baseName().split(sep), fileinfo.completeSuffix()]
			lls += [list([dir_path] + basename + [suffix])]
		lens = list(set([len(x) for x in lls]))
		if len(lens) != 1:
			print(f"{METHOD_NAME}: Input cannot be classified (different basename parts numbers).", end='\n')
			return False
		length = lens[0] - 2
		table = pd.DataFrame(lls, columns=['dir'] + list(range(lens[0] - 2)) + ['suffix'])
		if len(set(table['suffix'].to_list())) != 1:
			print(f"{METHOD_NAME}: Input cannot be classified (different suffixes).", end='\n')
			return False
		
		# --------------
		
		classes = pd.DataFrame(columns=['type', 'value'], index=range(length))
		
		val, typ = [[], []]
		for i in range(length):
			s = list(set(table[i].to_list()))
			val += ([s[0]] if len(s) == 1 else [s])
			typ += (['const'] if len(s) == 1 else [None])
		classes['value'], classes['type'] = [val, typ]
		
		
		
		print(table)
		print(classes)

Blister.Logo("Test!")

Huita.ClassifyInput(Blister.Input(['/dev/datasets/ngs_data/SeqDataTransfer/FatimannusaQadri/*.fastq.gz']))
