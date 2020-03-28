__version__ = "0.2"
__author__ = "regnveig"

import argparse
import os
import pandas
import pysam
import sys

METRICS_INDEX = ["Total, reads", "Unmapped, reads", "Strand Dups, reads", "Unmapped, %", "Strand Dups, %"]

def createParser():
	
	Default_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="Strandless: duplicate removal", epilog="Author: regnveig")
	Default_parser.add_argument('-v', '--version', action='version', version=__version__)
	Default_parser.add_argument ('-f', '--format', default='SAM', choices=['SAM', 'BAM'], dest="out_format", help='Output format (SAM, BAM). Default = SAM')
	Default_parser.add_argument ('-i', '--input', required=True, dest="in_stream", help='Input file. - for stdin')
	Default_parser.add_argument ('-o', '--output', default='-', dest="out_stream", help='Output file. Default = stdout')
	Default_parser.add_argument ('-m', '--metrics', default=None, dest="metrics_file", help='Metrics file. Default = None')
	
	return Default_parser

def pack_handle(pack):
	result = []
	sorter = {}
	for read in pack: 
		if read.is_unmapped: result += [read]
		else: sorter[str(read.reference_start) + ":" + str(read.reference_length)] = read
	return ((len(result), len(pack) - len(result) - len(sorter)), result + list(sorter.values()))

def strandless_func(namespace):
	
	_format = "w" if namespace.out_format == 'SAM' else "wb"
	metrics = pandas.Series([0, 0, 0, None, None], index=METRICS_INDEX)
	samfile = pysam.AlignmentFile(namespace.in_stream, "r")
	filtered = pysam.AlignmentFile(namespace.out_stream, _format, template=samfile)
	metrics_file = open(namespace.metrics_file, 'wt') if (namespace.metrics_file is not None) else None
		
	sorter = {}
	
	for read in samfile:
		metrics[METRICS_INDEX[0]] += 1
		try:
			sorter[read.query_name] += [read]
		except KeyError:
			if sorter:
				result = pack_handle(list(sorter.values())[0])
				metrics[METRICS_INDEX[1]] += result[0][0]
				metrics[METRICS_INDEX[2]] += result[0][1]
				for read2 in result[1]: filtered.write(read2)
				sorter = {}
			sorter[read.query_name] = [read]
			
	result = pack_handle(list(sorter.values())[0])
	metrics[METRICS_INDEX[1]] += result[0][0]
	metrics[METRICS_INDEX[2]] += result[0][1]
	for read2 in result[1]: filtered.write(read2)
	
	samfile.close()
	filtered.close()
	
	# report prepare
	metrics[METRICS_INDEX[3]] = metrics[METRICS_INDEX[1]] * 100 / metrics[METRICS_INDEX[0]]
	metrics[METRICS_INDEX[4]] = metrics[METRICS_INDEX[2]] * 100 / metrics[METRICS_INDEX[0]]
	report = "# Strandless Metrics\n# Input: " + namespace.in_stream + "\n# Output: " + namespace.out_stream + "\n" + metrics.to_frame().transpose().to_csv(sep='\t', index=False)
	if namespace.metrics_file is None: print(report, end='\n', file=sys.stderr)
	else:
		metrics_file.write(report)
		metrics_file.close()

if __name__ == '__main__':
	parser = createParser()
	namespace = parser.parse_args(sys.argv[1:])
	strandless_func(namespace)
