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

def blister_show():
	print('''
Blister v0.2

--------- I/O ---------

input_filenames = blister_input(filenames_list); if not input_filenames: exit()
output_dir = blister_dir(dir_path, create=True); if not output_dir: exit()
output_filename = blister_output(input_filename, output_dir, mod, suffix, rewrite=True); if not output_filename: continue
with blister_read(input_filename, mode='rt') as input_file:

--------- CLI ---------

blister_erase(n=1)
blister_progressbar(percent, start_time)
with blister_timestamp(title):

--------- CPU ---------

with blister_threading(title, THREADS_NUM = cpu_count()) as pool:
	results = pool.map(functools.partial(handler_func, var1=value1, var2=value2), iterable)
''')

def blister_input(filenames):
	file_list = []
	fileinfo_list = []
	fileinfo_unreadable = []
	for filename in filenames:
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
		print(f"Blister: List of unreadable files (will not be processed):", end="\n")
		for fileinfo in fileinfo_unreadable: print(f"\t{fileinfo}", end="\n")
	if not fileinfo_list:
		print(f"Blister: No input files exist.", end="\n")
		return False
	print(f"Blister: List of input files:", end="\n")
	for fileinfo in fileinfo_list: print(f"\t{fileinfo}", end="\n")
	return fileinfo_list

def blister_output(filename, output_dir, mod, suffix, rewrite=True):
	if mod != "": mod = "_" + mod
	if (suffix != ""): suffix = "." + suffix
	fileinfo_old = QFileInfo(filename)
	fileinfo = QFileInfo(QDir(output_dir), fileinfo_old.baseName() + mod + suffix)
	if (fileinfo.exists() and (not fileinfo.isFile())):
		print(f"Blister: This path is a dir:\n\t{fileinfo.absoluteFilePath()}", end="\n")
		return False
	if ((fileinfo.exists() and (not fileinfo.isWritable())) or ((not fileinfo.exists()) and (not QFileInfo(fileinfo.absolutePath()).permission(QFile.WriteUser)))):
		print(f"Blister: Writing this file is not permitted:\n\t{fileinfo.absoluteFilePath()}", end="\n")
		return False
	if (fileinfo.exists() and (rewrite == False)):
		fileinfo = QFileInfo(QDir(output_dir), fileinfo_old.baseName()+ "_" + str(int(time.time()) % 100000) + suffix)
		print(f"Blister: File to write already exists [rewriting is forbidden]. It will be renamed:\n\t{fileinfo_old.absoluteFilePath()} --> {fileinfo.absoluteFilePath()}", end="\n")
	return fileinfo.absoluteFilePath()

def blister_dir(_dir, create=True):
	dir_info = QFileInfo(_dir)
	if (dir_info.exists() and (not dir_info.permission(QFile.WriteUser))):
		print(f"Blister: Writing in this dir is not permitted:\n\t{dir_info.absoluteFilePath()}", end="\n")
		return False
	if ((not dir_info.exists()) and (not create)):
		print(f"Blister: This dir does not exist [creating new is forbidden]:\n\t{dir_info.absoluteFilePath()}", end="\n")
		return False
	if ((not dir_info.exists()) and create):
		result = QDir.mkpath(QDir(), dir_info.absoluteFilePath())
		if not result:
			print(f"Blister: Creating new dir was failed:\n\t{dir_info.absoluteFilePath()}", end="\n")
			return False
		else:
			print(f"Blister: New dir was created:\n\t{dir_info.absoluteFilePath()}", end="\n")
	return dir_info.absoluteFilePath()

def blister_gzip_check(filename):
	GZIP_MAGIC_NUMBER = "1f8b"
	with open(filename, 'rb') as file_check:
		return file_check.read(2).hex() == GZIP_MAGIC_NUMBER

def blister_bzip2_check(filename):
	BZIP2_MAGIC_NUMBER = "425a68"
	with open(filename, 'rb') as file_check:
		return file_check.read(3).hex() == BZIP2_MAGIC_NUMBER

@contextmanager
def blister_read(filename, mode='rt'):
	try:
		is_gz = blister_gzip_check(filename)
		is_bz2 = blister_bzip2_check(filename)
		file_format = "GZIP " if is_gz else ("BZIP2 " if is_bz2 else "")
		print(f"Blister: Open {file_format}file: {filename}", end="\n")
		f_obj = gzip.open(filename, mode) if is_gz else (bz2.open(filename, mode) if is_bz2 else open(filename, mode))
		yield f_obj
	except OSError:
		print(f"Blister: Unknown OSError: {filename}", end="\n")
	finally:
		f_obj.close()

def blister_sec2time(sec):
	return str(datetime.timedelta(seconds=int(sec)))

def blister_erase(n=1):	
	print("", end="\n")
	for _ in range(n):
		sys.stdout.write('\x1b[1A') # cursor up line
		sys.stdout.write('\x1b[2K') # erase line

def blister_progressbar(percent, start_time):
	fill = "#" * int(percent * 50)
	empty = "." * int(50 - len(fill))
	print(f"\t[{fill}{empty}] %4.0f%% [%s]" % (percent * 100, blister_sec2time(time.time() - start_time)), end="\r")

@contextmanager  
def blister_timestamp(title):
	start_time = time.time()
	yield
	print(f"Blister: {title} is done [%s]" % (blister_sec2time(time.time() - start_time)), end="\n")

@contextmanager
def blister_threading(title, THREADS_NUM = cpu_count()):
	print(f"Blister: starting {title} in {THREADS_NUM} cores ...", end="\n")
	pool = Pool(THREADS_NUM)
	yield pool
	pool.close()
	pool.join()
	del pool
	print(f"Blister: {title} threads are done.", end="\n")

'''
# TEMPLATES

## Pickle

import pickle

### dump
with open(dump_filename, 'wb') as f:
	pickle.dump(data, f)

### load
with open(dump_filename, 'rb') as f:
	data_new = pickle.load(f)

## Simple plot

import matplotlib.pyplot as plt

plt.plot(list)
plt.yscale("log"|"linear")
plt.ylabel('Y label')
plt.xlabel('X label')
plt.suptitle('Title')
plt.savefig("filename")

## PyVCF

import vcf

vcf_reader = vcf.Reader(open(vcf_filename, 'r'))

for record in vcf_reader:
	### do stuff

## BioPython

from Bio import SeqIO

### input
with open(input_filename, "rU") as input_handle:
	for record in SeqIO.parse(input_handle, input_format):
		# do stuff

### output
with open(output_filename, "w") as output_handle:
	SeqIO.write(sequences, output_handle, output_format)
	
### convert
count = SeqIO.convert(input_filename, input_format, output_filename, output_format)
print(f"Converted {count} records from {input_format} to {output_format}", end='\n')

## CSV

import pandas as pd

data = pd.read_csv(input_filename, sep=',', nrows=None)
data.to_csv(output_filename, sep='\t', index=False, mode='w')

'''
