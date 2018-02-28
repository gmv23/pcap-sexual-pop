#!/usr/bin/python

import sys

file1, file2 = sys.argv[1:]

fh = open(file1, "r")

file1_indvs = []

for line in fh:
	line = line.strip().split(":")
	file1_indvs.append(line[0])

fh.close()

fh = open(file2, "r")

for line in fh:
	line = line.strip().split(":")
	if line[0] not in file1_indvs:
		print line[0]

fh.close()
