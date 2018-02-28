#!/usr/bin/python

# The point of this script is to take samples that are matching in DNASample field in keyfile
# And give them the same FullSampleName field
# This field will now show each field from previous samples separated by a comma
# ie DNASample:Flowcell1,Flowcell2:Lane1,Lane2:LibraryPrepID1,LibraryPrepID2

import sys

input, output = sys.argv[1:]

sample_dic = {}

#Make dictionary of repeated samples

fh = open(input, "r")
fh.readline()

for line in fh:
	line = line.strip().split("\t")
	DNASample = line[3]
	extra_fields = [line[0], line[1], line[7]]
	if DNASample not in sample_dic.keys():
		sample_dic[DNASample] = []
		for field in extra_fields:
			sample_dic[DNASample].append([field])
	else:
		for i in range(0,len(sample_dic[DNASample])):
			sample_dic[DNASample][i].append(extra_fields[i])		

fh.close()

#Print infomration on how many times each sample is represented

counts_dic = {}

for sample in sample_dic.keys():
	counts = len(sample_dic[sample][0])
	if counts not in counts_dic.keys():
		counts_dic[counts] = 1
	else:
		counts_dic[counts] += 1

for counts in counts_dic.keys():
	print "%i samples represented %i times" %(counts_dic[counts], counts)

#Loop through key file again, make new fullsamplename for repeated samples and write new file

out = open(output, "w")
fh = open(input, "r")

out.write(fh.readline())

for line in fh:
	line = line.strip().split("\t")
	DNASample = line[3]
	if len(sample_dic[DNASample][1]) == 1:
		out.write("\t".join(line) + "\n")
	elif len(sample_dic[DNASample][1]) > 1:
		line[18] = DNASample + ":" + ":".join([",".join(field) for field in sample_dic[DNASample]])
		out.write("\t".join(line) + "\n")		
	else:
		print "problem..."



out.close()
fh.close()

