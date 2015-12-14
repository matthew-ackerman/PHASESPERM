#!/usr/bin/python

import sys
pairs={}
File=open(sys.argv[1])
for line in File:
	line=line.strip('\n').split(',')
	name=line[6]
	if name[0:5]=="GC(2)":
		MIN=int(line[3])
		MAX=int(line[4])
		name=name.split('#')
		SPERM=line[0]+"|"+name[0].split(')')[1]
		EVENT=name[1].split('_')[0]
		try:
			pairs[SPERM][EVENT][2]=MIN
			pairs[SPERM][EVENT][3]=MAX
		except:
			try:
				pairs[SPERM][EVENT]=[MIN, MAX, MIN, MAX]
			except:
				pairs[SPERM]={}
				pairs[SPERM][EVENT]=[MIN, MAX, MIN, MAX]
for s in pairs.keys():
	for e in pairs[s].keys():
		print pairs[s][e][2]-pairs[s][e][1]+1, '\b,' ,pairs[s][e][3]-pairs[s][e][0]+1
