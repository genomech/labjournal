
class Cuffcake(object):
	"""
My own library for Cufflinks results.
"""
	@staticmethod
	def ExpressionPerSample(cuff_dir, genes_id=None):
		"""
Return a table with genes expression per sample (Cuffdiff).
Can filter genes by list.
"""
		METHOD_NAME = f"Cuffcake.ExpressionPerSample"
		qdir = QDir(cuff_dir)
		tracking_fileinfo = QFileInfo(qdir, "genes.read_group_tracking")
		samples_fileinfo = QFileInfo(qdir, "read_groups.info")
		if (not tracking_fileinfo.exists()) or (not tracking_fileinfo.isFile()) or (not samples_fileinfo.exists()) or (not samples_fileinfo.isFile()):
			print(f"{METHOD_NAME}: Bad Cuffdiff directory:\n\t{qdir.absolutePath()}", end='\n')
			return False
		if (not tracking_fileinfo.isReadable()) or (not samples_fileinfo.isReadable()):
			print(f"{METHOD_NAME}: Cuffdiff files are not readable:\n\t{qdir.absolutePath()}", end='\n')
			return False
		tracking = pd.read_csv(tracking_fileinfo.absoluteFilePath(), sep='\t')
		samples = pd.read_csv(samples_fileinfo.absoluteFilePath(), sep='\t')
		samples_cols = samples[['file', 'condition', 'replicate_num']].applymap(str).apply(list, axis=1).to_list()
		samples_cols = [', '.join(x) for x in samples_cols]
		if genes_id != None:
			tracking = tracking[tracking['tracking_id'].isin(genes_id)]
			if tracking.shape[0] == 0:
				print(f"{METHOD_NAME}: No genes were found by list.", end='\n')
				return False
		samples.rename(columns={'replicate_num' : 'replicate'}, inplace=True)
		samples.drop(columns=['total_mass', 'norm_mass', 'internal_scale', 'external_scale'], inplace=True, axis=1)
		tracking = pd.merge(samples, tracking, on=['replicate', 'condition'], how='right')
		tracking.drop(columns=['raw_frags', 'internal_scaled_frags', 'external_scaled_frags', 'effective_length', 'status'], inplace=True, axis=1)
		tracking.sort_values(by=['tracking_id', 'condition', 'replicate'], ascending=True, inplace=True)
		cols = tracking.columns
		tracking = tracking.groupby('tracking_id').apply(pd.Series.to_list).apply(pd.DataFrame, columns=cols).apply(pd.DataFrame.drop, columns=['condition', 'replicate', 'tracking_id'], axis=1).to_dict()
		table = 0
		for it in tracking.keys():
			tracking[it] = tracking[it].rename(columns={'FPKM' : it})
			if type(table) == type(pd.DataFrame()): table = pd.merge(tracking[it], table, on=['file'], how='right')
			else: table = tracking[it]
		table = table.transpose()
		table.columns = samples_cols
		table.drop(labels='file', axis=0, inplace=True)
		return table

	@staticmethod
	def ReadMerged(filename):
		"""
Read merged.gtf (CuffMerge).
"""
		table = pd.read_csv(filename, sep='\t', header=None, names=['sequence', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
		header = [ (item.split(' ')[0]) for item in ((table['attributes'][0] + ' ').split('; ')[:-1]) ]
		attributes = []
		for attr in table['attributes'].items():
			attributes += [[ (item.split(' ')[1][1:-1]) for item in ((attr[1] + ' ').split('; ')[:-1]) ]]
		attributes = pd.DataFrame(attributes, columns=header)
		table = pd.concat([table, attributes], axis=1, join="outer")
		del attributes
		table.drop(columns=['attributes'], inplace=True)
		return table

	@staticmethod
	def LociIntersection(table, reference, table_format='one', reference_format='one'):
		"""
Returns a table of loci intersection.

_format:
1) one = {'LOCUS' : 'locus', 'CHROM' : None, 'POS' : None, 'BEGIN' : None, 'END' : None }
2) two = {'LOCUS' : None, 'CHROM' : 'chrom', 'POS' : 'pos', 'BEGIN' : None, 'END' : None }
3) three = {'LOCUS' : None, 'CHROM' : 'chrom', 'POS' : None, 'BEGIN' : 'begin', 'END' : 'end' }
"""
		METHOD_NAME = f"Cuffcake.LociIntersection"
		
		def pos2interval(pos):
			pos = pos.split('-')
			pos = [int(x) for x in pos]
			return I.IntInterval.closed(pos[0], pos[1])
		
		def prepare_table(_table, _format):
			_scheme = None
			if _format == 'one': _scheme = {'LOCUS' : 'locus', 'CHROM' : None, 'POS' : None, 'BEGIN' : None, 'END' : None }
			elif _format == 'two': _scheme = {'LOCUS' : None, 'CHROM' : 'chrom', 'POS' : 'pos', 'BEGIN' : None, 'END' : None }
			elif _format == 'three': _scheme = {'LOCUS' : None, 'CHROM' : 'chrom', 'POS' : None, 'BEGIN' : 'begin', 'END' : 'end' }
			else:
				if type(_format) != type(dict()):
					print(f"{METHOD_NAME}: Invalid format data.", end='\n')
					return False
				_scheme = _format
			
			if _scheme['LOCUS'] != None:
				if not (_scheme['LOCUS'] in _table):
					print(f"{METHOD_NAME}: Column '{_scheme['LOCUS']}' does not exist.", end='\n')
					return False
				# needs check string
				_table['123_chrom'] = _table[_scheme['LOCUS']].apply(lambda x: x.split(':')[0])
				_table['123_interval'] = _table[_scheme['LOCUS']].apply(lambda x: x.split(':')[1]).apply(pos2interval)
				_table.drop(columns=[_scheme['LOCUS']], inplace=True)
			else:
				if _scheme['CHROM'] != None:
					if not (_scheme['CHROM'] in _table):
						print(f"{METHOD_NAME}: Column '{_scheme['CHROM']}' does not exist.", end='\n')
						return False
					_table['123_chrom'] = _table[_scheme['CHROM']]
					_table.drop(columns=[_scheme['CHROM']], inplace=True)
					if _scheme['POS'] != None:
						if not (_scheme['POS'] in _table):
							print(f"{METHOD_NAME}: Column '{_scheme['POS']}' does not exist.", end='\n')
							return False
						_table['123_interval'] = _table[_scheme['POS']].apply(pos2interval)
						_table.drop(columns=[_scheme['POS']], inplace=True)
					elif (_scheme['BEGIN'] != None) and (_scheme['END'] != None):
						if not ((_scheme['BEGIN'] in _table) and (_scheme['END'] in _table)):
							print(f"{METHOD_NAME}: Column(s) '{_scheme['BEGIN']}', '{_scheme['END']}' do(es) not exist.", end='\n')
							return False
						_table['123_interval'] = _table[[_scheme['BEGIN'], _scheme['END']]].apply(lambda x: I.IntInterval.closed(int(x[_scheme['BEGIN']]), int(x[_scheme['END']])), axis=1)
						_table.drop(columns=[_scheme['BEGIN'], _scheme['END']], inplace=True)
					else:
						print(f"{METHOD_NAME}: Invalid format data.", end='\n')
						return False
				else:
					print(f"{METHOD_NAME}: Invalid format data.", end='\n')
					return False
			
			return _table

		table = prepare_table(table, table_format)
		if type(table) != type(pd.DataFrame()): return False
		reference = prepare_table(reference, reference_format)
		if type(reference) != type(pd.DataFrame()): return False

		def sort_func(a, b): return a['123_interval'].lower < b['123_interval'].lower
		def match_func(a, b): return not ((a['123_interval'].lower > b['123_interval'].upper) or (a['123_interval'].upper < b['123_interval'].lower))

		chroms = list(set(reference['123_chrom'].to_list()))
		sorted_genome = dict()

		for chrom in chroms:
			chr_table = reference[reference['123_chrom'] == chrom]
			sorted_genome[chrom] = BST(sort_func, match_func)
			for it in chr_table.iterrows(): sorted_genome[chrom].Insert(it[1])
			del chr_table
		
		ref_header = list(reference.columns.values)
		del reference
		
		def line_handle(row):
			lines = pd.DataFrame(sorted_genome[row[1]['123_chrom']].Search(row[1]), columns=ref_header) if (row[1]['123_chrom'] in sorted_genome.keys()) else pd.DataFrame(columns=ref_header, index=range(1))
			lines = lines if (lines.shape[0] != 0) else pd.DataFrame(columns=ref_header, index=range(1))
			lines['intersection'] = lines['123_interval'].apply(lambda x: (((row[1]['123_interval'] & x).upper - (row[1]['123_interval'] & x).lower) / (x.upper - x.lower)) if (x == x) else x)
			lines.rename(columns={'123_chrom': 'chrom'}, inplace=True)
			
			lines['start'] = lines['123_interval'].apply(lambda x: int(x.lower) if (x == x) else x)
			lines['end'] = lines['123_interval'].apply(lambda x: int(x.upper) if (x == x) else x)
			lines.drop(columns=['123_interval'], inplace=True)
			lines = lines.add_suffix("_ref")
			row = pd.DataFrame(([row[1]] * lines.shape[0]))
			row.rename(columns={'123_chrom': 'chrom'}, inplace=True)
			row['start'] = row['123_interval'].apply(lambda x: int(x.lower) if (x == x) else x)
			row['end'] = row['123_interval'].apply(lambda x: int(x.upper) if (x == x) else x)
			row.drop(columns=['123_interval'], inplace=True)
			row.index = lines.index.to_list()
			row = pd.concat([row, lines], axis=1)
			return row
		
		results = []
		for item in table.iterrows(): results += [line_handle(item)]
		big_table = pd.concat(results, ignore_index=True)
		
		return big_table
