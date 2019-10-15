#!/bin/python3

__author__ = "Anna Valeria [regnveig@yandex.ru]"
__date__ = "2019 Oct 7"
__version__ = "0.48"

from contextlib import contextmanager
from multiprocessing import cpu_count, Pool
from PyQt5.QtCore import QFileInfo, QFile, QDir
from tabulate import tabulate
import bz2
import datetime
import functools
import glob
import gzip
import intervals as I
import pandas as pd
import psutil
import subprocess
import sys
import time

class BST(object):
	"""
Primitive, but configurable binary search tree.
Useful if your class doesn't have search methods you need.

Usage:

def sort_func(a, b): return a > b
def match_func(a, b): return a == b
tree = Blister.BST(sort_func, match_func)
"""
	def __init__(self, sort_function, match_function):
		self.Sort_Function = sort_function
		self.Match_Function = match_function
		self.Root = None

	class Node: 
		def __init__(self, key): 
			self.left = None
			self.right = None
			self.val = key

	def Insert(self, value, root = -1):
		"""
Insert item to BST.

Usage: tree.Insert(value)
"""
		if root == -1: 
			node = self.Node(value)
			if self.Root is None: self.Root = node
			else:
				if self.Sort_Function(self.Root.val, node.val): 
					if self.Root.right is None: self.Root.right = node 
					else: self.Insert(value, self.Root.right) 
				else: 
					if self.Root.left is None: self.Root.left = node 
					else: self.Insert(value, self.Root.left)
		else:
			node = self.Node(value)
			if root is None: root = node
			else:
				if self.Sort_Function(root.val, node.val): 
					if root.right is None: root.right = node 
					else: self.Insert(value, root.right) 
				else: 
					if root.left is None: root.left = node 
					else: self.Insert(value, root.left)

	def Search(self, value, root = -1, lst = list()):
		"""
Return list of found items.

Usage: found = tree.Search(value)
"""
		if root == -1:
			if self.Root is not None:
				if self.Sort_Function(self.Root.val, value): lst = self.Search(value, self.Root.right, lst) 
				else: lst = self.Search(value, self.Root.left, lst)
				if self.Match_Function(self.Root.val, value): return lst + [self.Root.val]
			return lst
		else:
			if root is not None:
				if self.Sort_Function(root.val, value): lst = self.Search(value, root.right, lst) 
				else: lst = self.Search(value, root.left, lst)
				if self.Match_Function(root.val, value): return lst + [root.val]
			return lst

	def InOrder(self, root = -1):
		"""
Return sorted list of items.

Usage: sorted_list = tree.InOrder()
"""
		if root == -1:
			if self.Root:
				lst = self.InOrder(self.Root.left) 
				lst += [self.Root.val]
				lst += self.InOrder(self.Root.right)
				return lst
			else: return list()
		else:
			if root:
				lst = self.InOrder(root.left) 
				lst += [root.val]
				lst += self.InOrder(root.right)
				return lst
			else: return list()

class Blister(object):
	"""
Library for routine files and console procedures.
Based on KnV/GnP programming philosophy (Kostyl' & Velosiped, Govno & Palki).
Save my time for coffee, games, fap, and sleep.
"""
	@staticmethod
	def Input(filenames):
		"""
Take list of masks or filenames, return list of absolute filepaths, or False if error.
"""
		METHOD_NAME = f"Blister.Input"
		if type(filenames) != type(list()):
			print(f"{METHOD_NAME}: Invalid input type {type(filenames)}. List of strings only.", end='\n')
			return False
		file_list = []
		fileinfo_list = []
		fileinfo_unreadable = []
		for filename in filenames:
			if type(filename) != type(str()):
				print(f"{METHOD_NAME}: Invalid input type {type(filename)} in list. Strings only.", end='\n')
				return False
			file_list += glob.glob(QFileInfo(filename).absoluteFilePath())
		file_list = list(set(file_list))
		file_list.sort()
		for _file in file_list:
			fileinfo = QFileInfo(_file)
			if fileinfo.isFile():
				if fileinfo.isReadable():
					fileinfo_list += [fileinfo.absoluteFilePath()]
				else:
					fileinfo_unreadable += [fileinfo.absoluteFilePath()]
		if fileinfo_unreadable:
			print(f"{METHOD_NAME}: List of unreadable files (will not be processed):", end='\n')
			for fileinfo in fileinfo_unreadable: print(f"\t{fileinfo}", end='\n')
		if not fileinfo_list:
			print(f"{METHOD_NAME}: No input files exist or reachable.", end='\n')
			return False
		print(f"{METHOD_NAME}: List of input files:", end='\n')
		for fileinfo in fileinfo_list: print(f"\t{fileinfo}", end='\n')
		return fileinfo_list

	@staticmethod
	def Dir(dir_path, create=True):
		"""
Take dir path, return absolute path, or False if error.
Create new folder by default.
"""
		METHOD_NAME = f"Blister.Dir"
		if type(dir_path) != type(str()):
			print(f"{METHOD_NAME}: Invalid input type {type(dir_path)}. String only.", end='\n')
			return False
		dir_info = QFileInfo(dir_path)
		if (dir_info.exists() and (not dir_info.permission(QFile.WriteUser))):
			print(f"{METHOD_NAME}: Writing in this dir is not permitted:\n\t{dir_info.absoluteFilePath()}", end='\n')
			return False
		if ((not dir_info.exists()) and (not create)):
			print(f"{METHOD_NAME}: This dir does not exist [creating new is forbidden]:\n\t{dir_info.absoluteFilePath()}", end='\n')
			return False
		if ((not dir_info.exists()) and create):
			result = QDir.mkpath(QDir(), dir_info.absoluteFilePath())
			if not result:
				print(f"{METHOD_NAME}: Creating new dir was failed:\n\t{dir_info.absoluteFilePath()}", end='\n')
				return False
			else:
				print(f"{METHOD_NAME}: New dir was created:\n\t{dir_info.absoluteFilePath()}", end='\n')
		return dir_info.absoluteFilePath()

	@staticmethod
	def Output(filename, output_dir, mod, suffix, rewrite=True, index=-1):
		"""
Take input filename, output dir path, and make a new filename using mod and suffix.
Rewrite a file by default. Can be used with Blister.Threading().
"""
		METHOD_NAME = f"Blister.Output"
		thread_id = Blister.ThreadID(index)
		if mod != "": mod = "_" + mod
		if (suffix != ""): suffix = "." + suffix
		fileinfo_old = QFileInfo(filename)
		fileinfo = QFileInfo(QDir(output_dir), fileinfo_old.baseName() + mod + suffix)
		if (fileinfo.exists() and (not fileinfo.isFile())):
			print(f"{thread_id}{METHOD_NAME}: This path is a dir:\n{thread_id}\t{fileinfo.absoluteFilePath()}", end='\n')
			return False
		if ((fileinfo.exists() and (not fileinfo.isWritable())) or ((not fileinfo.exists()) and (not QFileInfo(fileinfo.absolutePath()).permission(QFile.WriteUser)))):
			print(f"{thread_id}{METHOD_NAME}: Writing this file is not permitted:\n{thread_id}\t{fileinfo.absoluteFilePath()}", end='\n')
			return False
		if (fileinfo.exists() and (rewrite == False)):
			fileinfo = QFileInfo(QDir(output_dir), fileinfo_old.baseName()+ "_" + str(int(time.time()) % 100000) + suffix)
			print(f"{thread_id}{METHOD_NAME}: File to write already exists [rewriting is forbidden]. It will be renamed:\n{thread_id}\t{fileinfo_old.absoluteFilePath()} --> {fileinfo.absoluteFilePath()}", end='\n')
		return fileinfo.absoluteFilePath()

	@staticmethod
	def GzipCheck(filename):
		"""
Check if file is gzipped.
"""
		GZIP_MAGIC_NUMBER = "1f8b"
		with open(filename, 'rb') as file_check:
			return file_check.read(2).hex() == GZIP_MAGIC_NUMBER

	@staticmethod
	def Bzip2Check(filename):
		"""
Check if file is bzipped.
"""
		BZIP2_MAGIC_NUMBER = "425a68"
		with open(filename, 'rb') as file_check:
			return file_check.read(3).hex() == BZIP2_MAGIC_NUMBER

	@staticmethod
	@contextmanager
	def Read(filename, mode='rt', index=-1):
		"""
Context manager for opening files. Can handle gzipped and bzipped too.
Yield opened file. Can be used with Blister.Threading().
"""
		METHOD_NAME = f"Blister.Read"
		try:
			thread_id = Blister.ThreadID(index)
			is_gz = Blister.GzipCheck(filename)
			is_bz2 = Blister.Bzip2Check(filename)
			f_obj = gzip.open(filename, mode) if is_gz else (bz2.open(filename, mode) if is_bz2 else open(filename, mode))
			yield f_obj
			f_obj.close()
		except OSError:
			print(f"{thread_id}{METHOD_NAME}: Unknown OSError:\n{thread_id}\t{filename}", end='\n')
			exit()

	@staticmethod
	def Seal(dir_path):
		"""
Seal a folder: count recursively MD5 checksum, write it to all.md5, then make the folder read-only.
Shell function, need package md5deep.
"""
		METHOD_NAME = f"Blister.Seal"
		qdir = QDir(dir_path)
		if not qdir.exists():
			print(f"{METHOD_NAME}: Path doesn't exist or isn't a dir.", end='\n')
			return False
		qfiles = QFileInfo(qdir, "**/*")
		file_list = glob.glob(qfiles.absoluteFilePath(), recursive=True)
		for _file in file_list:
			fileinfo = QFileInfo(_file)
			if not fileinfo.isReadable():
				print(f"{METHOD_NAME}: There is input cannot be checked (unreadable files). Stopped.", end='\n')
				return False
		start_time = time.time()
		sp = subprocess.Popen(f"cd {qdir.absolutePath()}; (md5deep -lr * > all.md5); chmod 555 -R {qdir.absolutePath()}", shell=True, stderr=subprocess.PIPE)
		out, err = sp.communicate()
		if err != b'':
			print(f"{METHOD_NAME}: Shell error: {str(err)}", end='\n')
			return False
		print(f"{METHOD_NAME}: Dir was sealed [%s]:\n\t{qdir.absolutePath()}" % (Blister.SecToTime(time.time() - start_time)), end='\n')
		return True

	@staticmethod
	def CountLines(filename):
		"""
Count lines if file. Can handle gzipped and bzipped too.
Shell function, need packages gzip and bzip2.
"""
		METHOD_NAME = f"Blister.CountLines"
		is_gz = Blister.GzipCheck(filename)
		is_bz2 = Blister.Bzip2Check(filename)
		command = f"zcat {filename} | wc -l" if is_gz else (f"bzip2 -dc {filename} | wc -l" if is_bz2 else f"cat {filename} | wc -l")
		sp = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = sp.communicate()
		if err != b'':
			print(f"{METHOD_NAME}: Shell error: {str(err)}", end='\n')
			return False
		return int(out)

	@staticmethod
	def Logo(script_name):
		"""
Just a beautiful minimalistic logo with timestamp.
Call Blister.Killer().
"""
		print(f"\n*** {script_name} ***\n\n{datetime.datetime.today().strftime('Now: %Y-%m-%d, %H:%M')}", end="\n\n")
		Blister.Killer()

	@staticmethod
	@contextmanager  
	def Timestamp(title, filename_1="", filename_2="", index=-1):
		"""
Context manager for timestamps. Can show input and output files if desired.
Yield start time. Can be used with Blister.Threading().
"""
		METHOD_NAME = f"Blister.Timestamp"
		thread_id = Blister.ThreadID(index)
		filename_1 = f":\n\t{filename_1}" if filename_1 != "" else f""
		filename_2 = f" --> {filename_2}" if filename_2 != "" else f""
		start_time = time.time()
		yield start_time
		print(f"{thread_id}{METHOD_NAME}: Process {title} is done [%s]{filename_1}{filename_2}" % (Blister.SecToTime(time.time() - start_time)), end='\n')
	
	@staticmethod
	def SecToTime(sec):
		"""
Return seconds in hh:mm:ss format.
"""
		return str(datetime.timedelta(seconds=int(sec)))

	@staticmethod
	def Erase(n=1):
		"""
Erase n last strings in console.
"""
		print("", end='\n')
		for _ in range(n):
			sys.stdout.write('\x1b[1A') # cursor up line
			sys.stdout.write('\x1b[2K') # erase line

	@staticmethod
	def ProgressBar(percent, start_time):
		"""
Draw a cute minimalistic progressbar.
Can take start time from Blister.Timestamp().
"""
		fill = "#" * int(percent * 50)
		empty = "." * int(50 - len(fill))
		print(f"\t[{fill}{empty}] %4.0f%% [%s]" % (percent * 100, Blister.SecToTime(time.time() - start_time)), end="\r")

	@staticmethod
	def Total(total, start_time):
		"""
Just a total number of something. Useless but cool.
Can take start time from Blister.Timestamp().
"""
		print(f"\tTotal: {total} [%s]" % (Blister.SecToTime(time.time() - start_time)), end="\r")

	@staticmethod
	@contextmanager
	def Threading(title, THREADS_NUM = cpu_count()):
		"""
Context manager for multiprocessing.

Usage:

def the_thread(block):
	index, item = block
	# do stuff
with Blister.Threading(title, THREADS_NUM = cpu_count()) as pool:
	results = pool.map(functools.partial(the_thread), enumerate(iterable))
"""
		METHOD_NAME = f"Blister.Threading"
		print(f"{METHOD_NAME}: Start {title} in {THREADS_NUM} cores ...", end='\n')
		start_time = time.time()
		pool = Pool(THREADS_NUM)
		yield pool
		pool.close()
		pool.join()
		del pool
		print(f"{METHOD_NAME}: {title} threads are done [%s]" % (Blister.SecToTime(time.time() - start_time)), end='\n')

	@staticmethod
	def EachFile(title, filenames, dir_path, THREADS_NUM = cpu_count()):
		"""
Decorator for threading by each file. Really cool stuff.

Usage:

def the_thread(block, output_dir):
	index, input_filename = block
	# do stuff
results = Blister.EachFile(title, [filenames], dir_path, THREADS_NUM = cpu_count())(the_thread)()
"""
		def Decorator(function):
			def Wrapper():
				Blister.Logo(title)
				input_filenames = Blister.Input(filenames)
				if not input_filenames: return False
				output_dir = Blister.Dir(dir_path, create=True)
				if not output_dir: return False
				with Blister.Threading(title, THREADS_NUM) as pool:
					results = pool.map(functools.partial(function, output_dir=output_dir), enumerate(input_filenames))
				return results
			return Wrapper
		return Decorator
	
	def ThreadID(index):
		return (f"[{index}]\t" if (index != -1) else f"")

	@staticmethod
	def Killer():
		"""
Kill all stopped processes and clean the memory.
"""
		METHOD_NAME = f"Blister.Killer"
		killed = pd.DataFrame(columns=['pid', 'name'])
		for proc in psutil.process_iter():
			try:
				pinfo = proc.as_dict(attrs=['pid', 'name', 'status'])
				if pinfo['status'] == 'stopped':
					proc.kill()
					killed = killed.append(pd.Series([pinfo['pid'], pinfo['name']], index=killed.columns), ignore_index=True)
			except psutil.NoSuchProcess:
				pass
		if killed.empty: return
		lst = set(killed['name'].to_list())
		print(f"{METHOD_NAME}: Stopped processes are killed:", end='\n')
		for victim in lst:
			print(f"\t{victim} - {killed[killed['name'] == victim]['pid'].to_list()}", end='\n')

	@staticmethod
	def Sleep(max_level=70.0, interval=1, index=-1):
		"""
Wait until memory percent is less than max level.
Can be useful if your script needs a lot of memory.
Can be used with Blister.Threading().
"""
		METHOD_NAME = f"Blister.Sleep"
		thread_id = Blister.ThreadID(index)
		if psutil.virtual_memory().percent <= max_level: return
		print(f"{thread_id}{METHOD_NAME}: Process is sleeping (memory = {psutil.virtual_memory().percent}%)", end='\n')
		while psutil.virtual_memory().percent > max_level:
			time.sleep(interval)

	@staticmethod
	def GitHubTable(dataframe, index=False):
		"""
My personal saviour. Turn pandas table into GitHub MD table.
Save lots of my time.
"""
		METHOD_NAME = f"Blister.GitHubTable"
		if type(dataframe) != type(pd.DataFrame()):
			print(f"{METHOD_NAME}: Invalid input type {type(dataframe)}. pandas.DataFrame only.", end='\n')
			return False
		return tabulate(dataframe, headers=dataframe.columns, tablefmt="github", showindex=index)

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

if __name__ == "__main__":
	sp = subprocess.Popen(f"pydoc3 -w ./Liebe.py; firefox ./Liebe.html", shell=True, stdout=subprocess.PIPE)
	out, err = sp.communicate()
