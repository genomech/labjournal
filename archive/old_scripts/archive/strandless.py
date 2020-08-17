#!/bin/python3

__version__ = "0.01"
__author__ = "Zanthia"

import argparse
import os
import sys

def createParser():
	Default_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="Strandless: duplicate removal", epilog="Author: regnveig")
	Default_parser.add_argument('-v', '--version', action='version', version=__version__)
	Default_parser.add_argument ('-f', '--format', default='SAM', choices=['SAM', 'BAM'], dest="out_format", help='Output format (SAM, BAM). Default = SAM')
	Default_parser.add_argument ('-i', '--input', required=True, dest="in_stream", help='Input file. - for stdin')
	Default_parser.add_argument ('-o', '--output', default='-', dest="out_stream", help='Output file. Default = stdout')
	Default_parser.add_argument ('-m', '--metrics', default=None, dest="metrics_file", help='Metrics file. Default = None')
	return Default_parser

def strandless_func(out_format, in_stream, out_stream, metrics_file):
	
	import pandas
	import pysam
	import random
	
	metrics = pandas.Series([0, 0, 0, None, None], index=["Total, reads", "Unmapped, reads", "Strand Dups, reads", "Unmapped, %", "Strand Dups, %"])
	samfile = pysam.AlignmentFile(in_stream, "r")
	filtered = pysam.AlignmentFile(out_stream, out_format, template=samfile)
	current_name = ""
	read_pack = []
	
	for read in samfile:
		metrics["Total, reads"] += 1
		if read.is_unmapped:
			filtered.write(read)
			metrics["Unmapped, reads"] += 1
			continue
		if read.query_name != current_name:
			for val in set([str(x.reference_start) + ":" + str(x.reference_length) for x in read_pack]):
				l = []
				for x in [(str(x.reference_start) + ":" + str(x.reference_length), x) for x in read_pack]:
					if x[0] == val: l += [x[1]]
				filtered.write(random.choice(l))
				metrics["Strand Dups, reads"] += len(l) - 1
			current_name = read.query_name
			read_pack = []
		read_pack += [read]
		
	metrics["Unmapped, %"] = metrics["Unmapped, reads"] * 100 / metrics["Total, reads"]
	metrics["Strand Dups, %"] = metrics["Strand Dups, reads"] * 100 / metrics["Total, reads"]
	with open(metrics_file, 'wt') as m:
		m.write("# Strandless Metrics\n# Input: " + in_stream + "\n# Output: " + out_stream + "\n" + metrics.to_frame().transpose().to_csv(sep='\t', index=False))
	samfile.close()
	filtered.close()

if __name__ == '__main__':
	parser = createParser()
	namespace = parser.parse_args(sys.argv[1:])
	_format = "wt" if (namespace.out_format == 'SAM') else "wb"
	strandless_func(_format, namespace.in_stream, namespace.out_stream, namespace.metrics_file)
