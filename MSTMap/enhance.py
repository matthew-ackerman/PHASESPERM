import sys
A=open(sys.argv[1])
B=open(sys.argv[2])
chm={}
for line in A:
	line=line.strip('\n').split()
	name=line[0].split('_')[1]
	loc=line[5]
	try:
		chm[name].append(loc)
	except:
		chm[name]=[loc]
for line in B:
	l=line.strip('\n').split()
	if len(l)==2:
		if l[0][0:2]=="sc":
			try:
				print '\t'.join(l)+'\t'+str(chm[l[0].split('_')[1] ])
			except:
				print '\t'.join(l)+"\t?"
		else:
			print '\t'.join(l)
	else:
		print '\t'.join(l)
