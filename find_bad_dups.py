#!/usr/bin/bash

#Print out file with lower depth duplicates

import sys

input, output = sys.argv[1:]

fh = open(input, "r")
fh.readline()

to_remove = []

depth_dic = {}

for line in fh:
	ind, depth = line.strip().split("\t")[0:3:2]
	ID = ind.strip().split(":")[0]
	if ID not in depth_dic.keys():
		depth_dic[ID] = [(ind, depth)]
	else:
		depth_dic[ID].append((ind,depth))

for key in depth_dic.keys():
	if len(depth_dic[key]) > 1:
		max_value = max([float(x[1]) for x in depth_dic[key]])
		for tup in depth_dic[key]:
			if float(tup[1]) < max_value:
				to_remove.append(tup[0])	

out = open(output, "w")
out.write("\n".join(to_remove))

fh.close()
