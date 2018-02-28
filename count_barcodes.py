#!/usr/bin/python

import os, sys

barcode_path, fastq_path, out_name = sys.argv[1:]

barcode_counts = {}
barcode_wells = {}
wells_in_order = []

fh = open(barcode_path, "r")
fh.readline() #pass header

for line in fh:
	barcode, well = line.strip().split("\t")
	barcode_counts[barcode] = 0
	barcode_wells[well] = barcode
	wells_in_order.append(well)

fh.close()

fh = open(fastq_path, "r")

line_count = 0
bad_read = 0

for line in fh:
	line_count += 1
	if line_count % 4 !=  2:
		continue
	seq = line.strip()
	cut_index = seq.find("CAGC")
	if cut_index == -1 or seq[:cut_index] not in barcode_counts.keys():
		bad_read += 1
	else:
		barcode_counts[seq[:cut_index]] = barcode_counts[seq[:cut_index]] + 1

fh.close()

out = open(out_name, "w")

out.write("bad read count:\t" + str(bad_read) + "\n")

for well in wells_in_order:
	out.write(well + "\t" + barcode_wells[well] + "\t" + str(barcode_counts[barcode_wells[well]]) + "\n")

out.close()
