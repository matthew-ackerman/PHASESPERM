#!/usr/bin/python
import sys
import random

N=['A','C','T','G']
Number_of_sperm=100
Number_of_scaffolds=100
REC=2
Linkagegroups=14

GC=20

Number_of_markers=50000
Prob_of_marker=0.50
Markers_per_scaffold=Number_of_markers/Number_of_scaffolds
Markers_per_scaffold=Number_of_markers/Number_of_scaffolds
PREC=float(REC*Linkagegroups)/float(Number_of_markers)
scaffolds=[]

File=open("test.txt", 'w')
micro=open("micro-test.txt", 'w')

for y in range(0, Linkagegroups):
	scaffolds.append([])
for x in range(0, Number_of_scaffolds):
	scaffolds[random.randint(0, Linkagegroups-1)].append("scaffold_"+str(x))
		
for y in range(0, Linkagegroups):
	phase=[]
	for s in range(0, Number_of_sperm):
		phase.append(random.randint(0, 1) )
	for scaf in scaffolds[y]:
		print scaf, " is in Linkagegroup ", y
		micro.write(scaf+"\t1\t2\tDp840\t100\t"+str(y)+"\t.\td155\n")
		for x in range(0, Markers_per_scaffold):
			line=[scaf, str(x*10)]
			n=N[:]
			P=[n.pop(random.randint(0, 3)), n.pop(random.randint(0, 2))]
			line.append(P[0])
			for s in range(0, Number_of_sperm):
				R=random.random()
				if(R<Prob_of_marker):
					line.append(P[phase[s]])
				else:
					line.append('-')
				R=random.random()
				if(R<PREC):
					phase[s]=int(phase[s]==0)
			File.write('\t'.join(line)+'\n')

