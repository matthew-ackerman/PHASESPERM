import sys
File=open(sys.argv[1])
x=3
for line in File:
	line=line.split()
	h=line[0:x]
	out=[line[x+6]]+line[x+8:]
	print '\t'.join(h+out)

