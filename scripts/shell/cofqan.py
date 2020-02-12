__version__ = "0.3"
__author__ = "regnveig"

from Bio import SeqIO
from contextlib import contextmanager
from multiprocessing import cpu_count, Pool
import argparse
import bz2
import functools
import gzip
import json
import re
import pandas
import string
import sys

GLOBAL_TABLE_COLS = ['Rate [%]', 'Mask']
GLOBAL_KMER_SIZE = 0
GLOBAL_UNRECOGNIZED = 'genome'
GLOBAL_MAX = 0
GLOBAL_THREADS = cpu_count()
GLOBAL_RATE_FLOOR = 0.01

def GzipCheck(filename):
	GZIP_MAGIC_NUMBER = "1f8b"
	with open(filename, 'rb') as file_check: return file_check.read(2).hex() == GZIP_MAGIC_NUMBER

def Bzip2Check(filename):
	BZIP2_MAGIC_NUMBER = "425a68"
	with open(filename, 'rb') as file_check: return file_check.read(3).hex() == BZIP2_MAGIC_NUMBER

@contextmanager
def Read(filename, mode='rt'):
	is_gz = GzipCheck(filename)
	is_bz2 = Bzip2Check(filename)
	f_obj = gzip.open(filename, mode) if is_gz else (bz2.open(filename, mode) if is_bz2 else open(filename, mode))
	yield f_obj
	f_obj.close()

@contextmanager
def Threading(threads):
	pool = Pool(threads)
	yield pool
	pool.close()
	pool.join()

def createParser():
	
	Default_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="CoFQan: Context FastQ Analyzer", epilog="Author: regnveig")
	Default_parser.add_argument('-v', '--version', action='version', version=__version__)
	Default_parser.add_argument ('-i', '--input', required=True, type=str, dest="input_filename", help='Fastq input file (may be gzipped or bzipped)')
	Default_parser.add_argument ('-o', '--output', required=True, type=str, dest="output_filename", help='Output CSV file')
	Default_parser.add_argument ('-p', '--patterns', required=True, type=str, dest="patterns", help='Patterns to look for, plain JSON format: \'{"pattern_1": "seq_1", "pattern_2": "seq_2"}\'')
	Default_parser.add_argument ('-k', '--kmer-size', default=GLOBAL_KMER_SIZE, type=int, dest="kmer_size", help=f'Max unrecognized K-mer size. Default = {GLOBAL_KMER_SIZE}')
	Default_parser.add_argument ('-u', '--unrecognized', default=GLOBAL_UNRECOGNIZED, type=str, dest="unrecognized", help=f'Long unrecognized sequences replacement. Default = {GLOBAL_UNRECOGNIZED}')
	Default_parser.add_argument ('-m', '--max-reads', default=GLOBAL_MAX, type=int, dest="max_reads", help=f'Max reads number to analyze (0 - no limit). Default = {GLOBAL_MAX}')
	Default_parser.add_argument ('-f', '--rate-floor', default=GLOBAL_RATE_FLOOR, type=float, dest="rate_floor", help=f'Min rate to write pattern in the table. Default = {GLOBAL_RATE_FLOOR}')
	Default_parser.add_argument ('-@', '--threads', default=GLOBAL_THREADS, type=int, dest="threads", help=f'Threads number. Default = {GLOBAL_THREADS}')
	
	return Default_parser

def seq_prepare(namespace, string):
	
	string = string.upper()
	for item in namespace.patterns.items(): string = string.replace(item[1], f"-{item[0]}-")
	for it in range(1, namespace.kmer_size): string = re.sub("(^|-)[ATGCN]{%s}($|-)" % (it), "--[%s]--" % it, string)
	
	return re.sub(r"[ATGCN]+", f"-{namespace.unrecognized}-", string)

def seq_count(main_list, string): return pandas.Series([main_list.count(string), string], index=GLOBAL_TABLE_COLS)

def read_context(namespace):
	
	main_list = []
	total = 0

	with Read(namespace.input_filename) as input_file:
		for record in SeqIO.parse(input_file, "fastq"):
			main_list += [record.seq.__str__()]
			total += 1
			if total == namespace.max_reads: break
	
	with Threading(namespace.threads) as pool1: main_list = pool1.map(functools.partial(seq_prepare, namespace), main_list)
	
	shorted_list = list(set(main_list))
	
	with Threading(namespace.threads) as pool2: main_table = pool2.map(functools.partial(seq_count, main_list), shorted_list)
		
	main_table = pandas.DataFrame(main_table)
	main_table[GLOBAL_TABLE_COLS[0]] = main_table[GLOBAL_TABLE_COLS[0]].apply(lambda x: x / (total / 100))
	main_table[main_table[GLOBAL_TABLE_COLS[0]] >= namespace.rate_floor].sort_values(by=[GLOBAL_TABLE_COLS[0]], ascending=False).to_csv(namespace.output_filename, sep=',', index=False)

if __name__ == '__main__':
	
	parser = createParser()
	namespace = parser.parse_args(sys.argv[1:])
	
	try:
		namespace.patterns = json.loads(namespace.patterns)
	except:
		print(f"Error: Patterns are not JSON format", end='\n', file=sys.stderr)
		exit(1)
	
	read_context(namespace)
