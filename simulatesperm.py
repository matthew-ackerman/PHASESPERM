#!/usr/bin/python
import sys
import random
import numpy

N=['A','C','T','G']						# Possible bases
Number_of_sperm=100						# Nuber of sperm sequenced
Number_of_scaffolds=100						# Number of scaffold in the genome
REC=2								# ?
Linkagegroups=14						# Number of linkage groups in the genome
heterozygosity=0.01						# Fraction of sites that are heterozygous

Number_of_markers=50000						# Number of markers in the file
Prob_of_marker=0.50						# The probabitly of a marker
Prob_of_error=0.0						# The probability that a single sperm is miscalled.

Prob_of_bad_start=0.0						# The probability of entering a region within a scaffold that is ?
Prob_of_bad_stop=0.1						# The probability of exiting a region within a scaffold that is ?

Prob_of_gc_start=0.9						# The probability of entering a gc region.

Markers_per_scaffold=Number_of_markers/Number_of_scaffolds	#
Markers_per_scaffold=Number_of_markers/Number_of_scaffolds	#
PREC=float(REC*Linkagegroups)/float(Number_of_markers)		#
scaffolds=[]							#

File=open("test.txt", 'w')
micro=open("micro-test.txt", 'w')
gctracks=open("gcpairs-test.txt", 'w')

isbad=False
isgctrack=False
writegc=False
gc_extention_probability=0.9

for y in range(0, Linkagegroups):
	scaffolds.append([])
for x in range(0, Number_of_scaffolds):
	scaffolds[random.randint(0, Linkagegroups-1)].append("scaffold_"+str(x))
		
for y in range(0, Linkagegroups):
	phase=[]
	for s in range(0, Number_of_sperm):
		phase.append(random.randint(0, 1) )
	for scaf in scaffolds[y]:
		pos=0
		print scaf, " is in Linkagegroup ", y
		micro.write(scaf+"\t1\t2\tDp840\t100\t"+str(y)+"\t.\td155\n")
		for x in range(0, Markers_per_scaffold):
			inter_marker_distance=numpy.random.geometric(heterozygosity)
			temp_marker_distance=inter_marker_distance
			pos+=inter_marker_distance
			line=[scaf, str(pos)]
			n=N[:]
			P=[n.pop(random.randint(0, 3)), n.pop(random.randint(0, 2))]
			line.append(P[0])
			if (writegc):
				maxgc+=inter_marker_distance
				gctracks.write(str(mingc+1)+", "+str(maxgc)+", "+str(truegc)+'\n' )
				writegc=False
				print "Writing GC"
			if not(isbad):
				if(random.random()<Prob_of_bad_start):
					badphase=[]
					for b in range(0, Number_of_sperm):
						badphase.append(random.randint(0, 1) )
						isbad=True
			if not(isgctrack):
				if(random.random()<Prob_of_gc_start):
					print "Starting GC"
					gcsperm=random.randint(0, Number_of_sperm-1)
					phase[gcsperm]=int(phase[gcsperm]==0)
					isgctrack=True
					writegc=False
					temp_marker_distance=random.randint(0, inter_marker_distance)
					maxgc=0
					mingc=0
					truegc=0
			if(isgctrack):
				Prob_of_gc_stop=1-gc_extention_probability**(temp_marker_distance)
				maxgc+=inter_marker_distance
				if(random.random()<Prob_of_gc_stop):
					truegc+=random.randint(0, temp_marker_distance)
					isgctrack=False
					print "Closing GC"
					phase[gcsperm]=int(phase[gcsperm]==0)
					writegc=True
				else:
					truegc+=temp_marker_distance
					mingc+=inter_marker_distance
			for s in range(0, Number_of_sperm):
				if (isbad):
					if(random.random()<Prob_of_marker):
						line.append(P[badphase[s]])
					else:
						line.append('-')
				else:
					if(random.random()<Prob_of_marker):
						if(random.random()>Prob_of_error):
							line.append(P[phase[s]])
						else:
							line.append(P[random.randint(0, 1)])
					else:
						line.append('-')
				R=random.random()
				if(R<PREC):
					phase[s]=int(phase[s]==0)
			File.write('\t'.join(line)+'\n')
			if(isbad):
				if(random.random()<Prob_of_bad_stop):
					isbad=False
