#!/usr/bin/python

#This script takes a directory of key files and combines them, only taking capsici or blank wells

import os, sys

path,outname = sys.argv[1:]

systemcall1 = "ls %s > key_files.txt" %path
os.system(systemcall1)

fh = open("./key_files.txt", "r")

files = []
for line in fh:
	files.append(line.strip())

fh.close()
os.system("rm key_files.txt")

file_count = 0
out = open(outname, "w")

for file in files:
	file_count += 1
	file_path = path + file
	fh = open(file_path, "r")
	if file_count != 1:
		fh.readline()
	for line in fh:
		species = line.strip().split("\t")[14]
		if species.upper() in [x.upper() for x in ["capsici", "na", "species", "blank", ""]]:
			out.write(line)
	fh.close()

out.close()

