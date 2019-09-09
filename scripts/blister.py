# Blister v0.3

from contextlib import contextmanager
from multiprocessing import cpu_count, Pool
from PyQt5.QtCore import QFileInfo, QFile, QDir
import bz2
import datetime
import functools
import glob
import gzip
import sys
import time

# --------- ETC ---------

# tabulate
# pandas: read_csv, to_csv
# Bio, vcf
# subprocess
# pickle
# matplotlib.pyplot

# def blister_seal(mask):

# --------- I/O ---------

# input_filenames = blister_input([filenames])
# if not input_filenames: exit()

def blister_input(filenames):
	if type(filenames) != type(list()):
		print(f"Blister: Invalid input type {type(filenames)}. List of strings only.", end='\n')
		return False
	file_list = []
	fileinfo_list = []
	fileinfo_unreadable = []
	for filename in filenames:
		if type(filename) != type(str()):
			print(f"Blister: Invalid input type {type(filename)} in list. Strings only.", end='\n')
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
		print(f"Blister: List of unreadable files (will not be processed):", end='\n')
		for fileinfo in fileinfo_unreadable: print(f"\t{fileinfo}", end='\n')
	if not fileinfo_list:
		print(f"Blister: No input files exist or reachable.", end='\n')
		return False
	print(f"Blister: List of input files:", end='\n')
	for fileinfo in fileinfo_list: print(f"\t{fileinfo}", end='\n')
	return fileinfo_list

# output_dir = blister_dir(dir_path, create=True)
# if not output_dir: exit()

def blister_dir(dir_path, create=True):
	if type(dir_path) != type(str()):
		print(f"Blister: Invalid input type {type(dir_path)}. String only.", end='\n')
		return False
	dir_info = QFileInfo(dir_path)
	if (dir_info.exists() and (not dir_info.permission(QFile.WriteUser))):
		print(f"Blister: Writing in this dir is not permitted:\n\t{dir_info.absoluteFilePath()}", end='\n')
		return False
	if ((not dir_info.exists()) and (not create)):
		print(f"Blister: This dir does not exist [creating new is forbidden]:\n\t{dir_info.absoluteFilePath()}", end='\n')
		return False
	if ((not dir_info.exists()) and create):
		result = QDir.mkpath(QDir(), dir_info.absoluteFilePath())
		if not result:
			print(f"Blister: Creating new dir was failed:\n\t{dir_info.absoluteFilePath()}", end='\n')
			return False
		else:
			print(f"Blister: New dir was created:\n\t{dir_info.absoluteFilePath()}", end='\n')
	print(f"Blister: Output dir has been chosen:\n\t{dir_info.absoluteFilePath()}", end='\n')
	return dir_info.absoluteFilePath()

# output_filename = blister_output(input_filename, output_dir, mod, suffix, rewrite=True)
# if not output_filename: continue

def blister_output(filename, output_dir, mod, suffix, rewrite=True, index=-1):
	thread_id = blister_thread_id(index)
	if mod != "": mod = "_" + mod
	if (suffix != ""): suffix = "." + suffix
	fileinfo_old = QFileInfo(filename)
	fileinfo = QFileInfo(QDir(output_dir), fileinfo_old.baseName() + mod + suffix)
	if (fileinfo.exists() and (not fileinfo.isFile())):
		print(f"{thread_id}Blister: This path is a dir:\n{thread_id}\t{fileinfo.absoluteFilePath()}", end='\n')
		return False
	if ((fileinfo.exists() and (not fileinfo.isWritable())) or ((not fileinfo.exists()) and (not QFileInfo(fileinfo.absolutePath()).permission(QFile.WriteUser)))):
		print(f"{thread_id}Blister: Writing this file is not permitted:\n{thread_id}\t{fileinfo.absoluteFilePath()}", end='\n')
		return False
	if (fileinfo.exists() and (rewrite == False)):
		fileinfo = QFileInfo(QDir(output_dir), fileinfo_old.baseName()+ "_" + str(int(time.time()) % 100000) + suffix)
		print(f"{thread_id}Blister: File to write already exists [rewriting is forbidden]. It will be renamed:\n{thread_id}\t{fileinfo_old.absoluteFilePath()} --> {fileinfo.absoluteFilePath()}", end='\n')
	print(f"{thread_id}Blister: File to write has been chosen:\n{thread_id}\t{fileinfo_old.absoluteFilePath()} --> {fileinfo.absoluteFilePath()}", end='\n')
	return fileinfo.absoluteFilePath()

# is_gzip = blister_gzip_check(filename)

def blister_gzip_check(filename):
	GZIP_MAGIC_NUMBER = "1f8b"
	with open(filename, 'rb') as file_check:
		return file_check.read(2).hex() == GZIP_MAGIC_NUMBER

# is_bzip2 = blister_bzip2_check(filename)

def blister_bzip2_check(filename):
	BZIP2_MAGIC_NUMBER = "425a68"
	with open(filename, 'rb') as file_check:
		return file_check.read(3).hex() == BZIP2_MAGIC_NUMBER

# with blister_read(input_filename, mode='rt') as input_file:

@contextmanager
def blister_read(filename, mode='rt', index=-1):
	try:
		thread_id = blister_thread_id(index)
		is_gz = blister_gzip_check(filename)
		is_bz2 = blister_bzip2_check(filename)
		file_format = "GZIP " if is_gz else ("BZIP2 " if is_bz2 else "")
		print(f"{thread_id}Blister: Open {file_format}file:\n{thread_id}\t{filename}", end='\n')
		f_obj = gzip.open(filename, mode) if is_gz else (bz2.open(filename, mode) if is_bz2 else open(filename, mode))
		yield f_obj
		f_obj.close()
	except OSError:
		print(f"{thread_id}Blister: Unknown OSError:\n{thread_id}\t{filename}", end='\n')
		exit()

# --------- CLI ---------

# blister_logo(script_name)

def blister_logo(script_name):
	print(f"\n*** {script_name} ***\n\n{datetime.datetime.today().strftime('Now: %Y-%m-%d, %H:%M')}", end="\n\n")

# with blister_timestamp(title, index) as start_time:

@contextmanager  
def blister_timestamp(title, index=-1):
	thread_id = blister_thread_id(index)
	start_time = time.time()
	yield start_time
	print(f"{thread_id}Blister: Process {title} is done [%s]" % (blister_sec2time(time.time() - start_time)), end='\n')

# blister_sec2time(sec)

def blister_sec2time(sec):
	return str(datetime.timedelta(seconds=int(sec)))

# blister_erase()

def blister_erase(n=1):	
	print("", end='\n')
	for _ in range(n):
		sys.stdout.write('\x1b[1A') # cursor up line
		sys.stdout.write('\x1b[2K') # erase line

# blister_progressbar(percent, start_time)

def blister_progressbar(percent, start_time):
	fill = "#" * int(percent * 50)
	empty = "." * int(50 - len(fill))
	print(f"\t[{fill}{empty}] %4.0f%% [%s]" % (percent * 100, blister_sec2time(time.time() - start_time)), end="\r")

# blister_total(total, start_time)

def blister_total(total, start_time):
	print(f"\tTotal: {total} [%s]" % (blister_sec2time(time.time() - start_time)), end="\r")

# --------- CPU ---------

# with blister_threading(title) as pool:
	# results = pool.map(functools.partial(the_thread), enumerate(iterable))
# def the_thread(block):
	# index = block[0]
	# item = block[1]

@contextmanager
def blister_threading(title, THREADS_NUM = cpu_count()):
	print(f"Blister: Start {title} in {THREADS_NUM} cores ...", end='\n')
	start_time = time.time()
	pool = Pool(THREADS_NUM)
	yield pool
	pool.close()
	pool.join()
	del pool
	print(f"Blister: {title} threads are done [%s]" % (blister_sec2time(time.time() - start_time)), end='\n')

# thread_id = blister_thread_id(index)

def blister_thread_id(index):
	return (f"[{index}]\t" if (index != -1) else f"")
