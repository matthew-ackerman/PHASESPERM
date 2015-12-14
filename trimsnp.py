#!/usr/bin/python
import sys

try:
	File=open(sys.argv[1])
	T=int(sys.argv[2])
except:
	print "Usage: "+sys.argv[0]+" <file.snp> <min calls>"
	quit()

for line in File:
	line=line.split()
	calls=''.join(line[3:])
	C=[calls.count('A'), calls.count('C'), calls.count('G'), calls.count('T')]
	C=sorted(C)
	if C[2]>T and C[3]>T:
		print '\t'.join(line)
