#!/usr/bin/bash

#Print out file with dups from key file

import sys

input = sys.argv[1]

sample_dic = {}

fh = open(input, "r")
fh.readline()

for line in fh:
	line = line.strip().split("\t")
	flowcell, sample, prepID, plateID, barcode = [line[0], line[3], line[7], line[8], line[2]]
	if sample not in sample_dic.keys():	
		sample_dic[sample] = [(flowcell, prepID, plateID, barcode)]
	else:
		sample_dic[sample].append((flowcell, prepID, plateID, barcode))

fh.close()

for sample in sample_dic.keys():
	if len(sample_dic[sample]) > 1:
		for cell in sample_dic[sample]:
			print sample + "\t" + "\t".join(cell)

