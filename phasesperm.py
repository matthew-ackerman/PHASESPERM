#!/usr/bin/python

import sys
import os
import numpy
from scipy.stats import poisson
import math
import operator

max_rec=15	
group_thresh=20
snp_thresh=5
minphased=90
N50=101
CGSIZE=10000
keep=False

try:
	for x in range(1, len(sys.argv)-1, 2):
		if sys.argv[x]=="-E":
			max_rec=int(sys.argv[x+1])
		if sys.argv[x]=="-G":
			group_thresh=int(sys.argv[x+1])
		if sys.argv[x]=="-S":
			snp_thresh=int(sys.argv[x+1])
		if sys.argv[x]=="-M":
			minphased=int(sys.argv[x+1])
		if sys.argv[x]=="-N50":	
			N50=int(sys.argv[x+1])
		if sys.argv[x]=="-h":
			raise Exception("help!")
		if sys.argv[x]=="-k":
			keep=True
	filename=sys.argv[-1]
	print "phasing "+filename+". This may take awhile."
	SpermFile=open(filename).read().split('\n')

except:
	print
	print "Usage: "+sys.argv[0]+" [OPTIONS] FILENAME"
	print
	print "NOTE: If the program generates an error of the form:"
	print "Traceback (most recent call last):"
	print "  File \"./SpermPhase-final3.py\", line 384, in <module>"
	print "    if int(this.split('_')[1] )<N50:"
	print "IndexError: list index out of range"
	print 
	print "try running it with the -N50 0 option."
	print
	print "\t-G : markers with fewer recombination events than this threshold are assumed to belong to the same linkage group and are phased accordingly."
	print "\t\tIf this threshold is set two low then markers can be miss-phased due to chance associations between linkage groups, and if it is set to high then markers"
	print "\t\tare not phased properly for the mapping program. Generally a value of 15% of your total sample works well, depending on the density of your markers."
	print
	print "\t-E : The maximum number of recombination events that can occur between adjesent marker/haplo-block in the reference genome. If a marker exceeds this"
	print "\t\t it is thrown out. Setting the number too low can result in throwing at all of your markers, and setting the number too high can result in including many "
	print "\t\t markers that are incorrectly phased because of some kind of error. Generally a value of around 15% will prevent poorly phased markers from being included in "
	print "\t\t the output, but allow most scaffolds to be phased."
	print
	print "\t-S : The minimum number of SNPs within a marker/haplo-block neccisary to attempt to use the marker in genetic map construction. The value of this number" 
	print "\t\tshould be set to a value large enough to exclude small markers/haplo-blocks that can arise from sequencing and alignment errors, but to a value small enough to"
	print "\t\trecover at least 2 x Sperm x Linkage Group markers."
	print 
	print "\t-M : The minimum number of sperm phased within a marker/haplo-block. Markers/haplo-block with less than this number of phased sperm are not used in any analysis."
	print
	print "\t-N50 : An integer reperesnting the largest number scaffold you wish to analalize. This generally should represent the portion of the genome that is relatively free "
	print "\t\t of assembly errors and repetive sequecnce. To disable this option entirely use 0. "
	print 
	print "\t-k : (keep) a flag to indicate that low coverage samples should be used in the analysis."
	print  
	print  
	print "Two output files will be generated : outR.csv which will print information about each SNP in the input file, and MST-data.txt which can be used by the program MSTMap"
	print "http://alumni.cs.ucr.edu/~yonghui/mstmap.html to produce your genetic map. Markers in MSTMap are named by ..."
	print 
	print "To finsh construction of the genetic map, complie MSTMap and type"
	print "MSTMap"
	print "This script depends on numpy and scipy. You will need to install these by..."
	print "on many linux distributions you can obtain these by typing"
	print "sudo apt-get numpy"
	print "and"
	print "sudo apt-get scipy"

	quit()
	
class recombination:
	def __init__ (self, line, start, stop, fivepSNPs, count):
		self.line=line
		self.start=start
		self.stop=stop
		self.end=stop
		self.count=1#count
		self.event=0
		self.fivepSNPs=fivepSNPs
		self.threepSNPs=1
		self.bad=False

def overlap (event1, event2):
	X=0
	if event1.fivepSNPs>X and event1.threepSNPs>X and event2.fivepSNPs>X and event2.threepSNPs>X:
		if (event2.start>event1.start):	
			if(event2.start<event1.stop):
				return True
		elif (event2.stop>event1.start):
			if(event2.stop<event1.stop):
				return True
		if (event1.start>event2.start):	
			if(event1.start<event2.stop):
				return True
		elif (event1.stop>event2.start):
			if(event1.stop<event2.stop):
				return True
	return False

def error (event, numlines, rate, P):
	if event.count>3:
		return True
	return False
	if ( (1.0-poisson.cdf(event.count-1, float(rate*numlines)*(event.stop-event.start) ) )<P ):
		return True
	return False

def bitCount(int_type):
	count = 0
	while(int_type):
		int_type &= int_type - 1
		count += 1
	return(count)

class SNP:
	def __init__(self, bp, calls, setbits):
		S=0
		M=0
		A='N'
		B='N'
		for x in xrange(0, len(calls) ):
			if A=='N':
				if (calls[x]!='-'):
					A=calls[x]
			elif B=='N':
				if (calls[x]!='-'):
					if calls[x]!=A:
						B=calls[x]
			M=M | ( (calls[x]=='-')<<x )
			S=S | ( (calls[x]==A)<<x )
		self.A=A
		self.B=B
		self.bp=bp
		self.M=M
		self.S=S
		self.setbits=setbits
		
	def flip(self):
		self.S=self.setbits&(~self.S)
		self.A, self.B = self.B, self.A
		return self
		

def getR(that, this):
	R1=this.setbits&( (~that.M)&( (~this.M)&(that.S^ this.S) ) )
	R2=this.setbits&( (~that.M)&( (~this.M)&(that.S^~this.S) ) )
	return [bitCount(R1), bitCount(R2), R1, R2]

class SNP_BLOCK:
	def __init__(self, setbits, numlines, scaffold):
		self.list=[]
		self.subblock=[]
		self.M=setbits
		self.S=0
		self.setbits=setbits
		self.numlines=numlines
		self.index={}
		self.scaffold=scaffold
		self.masked=False
		self.count=[]
		for x in range(0, numlines):
			self.count.append(0)
		
	def add_snp(self, this):
		[cR1, cR2, R1, R2]=getR(self, this)
		if cR1==0 and cR2==0 and len(self.list)>0:
			self.subblock.append(this)
			return False
		if cR1<cR2:
			if (cR1==0) :
				self.M=self.setbits&( ~( (~self.M)|(~this.M) ) )
				self.S=self.setbits&( (self.S&this.M) | (~this.M & this.S) )
				self.list.append(this)
				for x in range(0, numlines):
					self.count[x]+=(~this.M & 2**x!=0)
				for block in self.subblock: 
					[cR1, cR2, R1, R2]=getR(self, block)
					if cR1==0 and cR2>0:
						self.M=self.setbits&( ~( (~self.M)|(~block.M) ) )
						self.S=self.setbits&( (self.S&block.M) | (~block.M & block.S) )
						for x in range(0, numlines):
							self.count[x]+=(~block.M & 2**x!=0)
						self.list.append(block)
						self.subblock.remove(block)
					elif cR2==0 and cR1>0:
						block.flip()
						self.M=self.setbits&( ~( (~block.M)|(~block.M) ) )
						self.S=setbits&( (self.S&block.M) | (~block.M & block.S) )
						for x in range(0, numlines):
							self.count[x]+=(~block.M & 2**x!=0)
						self.list.append(block)
						self.subblock.remove(block)
					else:
						return True
				return False
			else:
				for block in self.subblock:
					[cR1, cR2, R1, R2]=getR(self, block)
					if cR1==0 and cR2>0:
						self.M=self.setbits&( ~( (~self.M)|(~block.M) ) )
						self.S=self.setbits&( (self.S&block.M) | (~block.M & block.S) )
						for x in range(0, numlines):
							self.count[x]+=(~block.M & 2**x!=0)
						self.list.append(block)
					elif cR2==0 and cR1>0:
						block.flip()
						self.M=self.setbits&( ~( (~block.M)|(~block.M) ) )
						self.S=setbits&( (self.S&block.M) | (~block.M & block.S) )
						for x in range(0, numlines):
							self.count[x]+=(~block.M & 2**x!=0)
						self.list.append(block)     
				return True
		else:
			if (cR2==0) :
				this.flip()
				self.M=self.setbits&( ~( (~self.M)|(~this.M) ) )
				self.S=setbits&( (self.S&this.M) | (~this.M & this.S) )
				for x in range(0, numlines):
					self.count[x]+=(~this.M & 2**x!=0)
				self.list.append(this)
				for block in self.subblock: 
					[cR1, cR2, R1, R2]=getR(self, block)
					if cR1==0 and cR2>0:
						self.M=self.setbits&( ~( (~self.M)|(~block.M) ) )
						self.S=self.setbits&( (self.S&block.M) | (~block.M & block.S) )
						for x in range(0, numlines):
							self.count[x]+=(~block.M & 2**x!=0)
						self.list.append(block)
						self.subblock.remove(block)
					elif cR2==0 and cR1>0:
						block.flip()
						self.M=self.setbits&( ~( (~block.M)|(~block.M) ) )
						self.S=setbits&( (self.S&block.M) | (~block.M & block.S) )
						for x in range(0, numlines):
							self.count[x]+=(~block.M & 2**x!=0)
						self.list.append(block)
						self.subblock.remove(block)
					else:
						return True
				return False
			else:
				for block in self.subblock:
					[cR1, cR2, R1, R2]=getR(self, block)
					if cR1==0 and cR2>0:
						self.M=self.setbits&( ~( (~self.M)|(~block.M) ) )
						self.S=self.setbits&( (self.S&block.M) | (~block.M & block.S) )
						for x in range(0, numlines):
							self.count[x]+=(~block.M & 2**x!=0)
						self.list.append(block)
					elif cR2==0 and cR1>0:
						block.flip()
						self.M=self.setbits&( ~( (~block.M)|(~block.M) ) )
						self.S=setbits&( (self.S&block.M) | (~block.M & block.S) )
						for x in range(0, numlines):
							self.count[x]+=(~block.M & 2**x!=0)
						self.list.append(block)
				this.flip()
				return True
	def flip(self):
		self.S=self.setbits&(~self.S)
		newlist=[]
		for snp in self.list:
			newlist.append(snp.flip() )
		self.list=newlist
		return self

	def write(self, scaffold, name, bR1, bR2):
		ret=[]
		for snp in self.list:
			bit=[]
			for x in range(0, self.numlines):
				if (~(snp.M)&(2**x))!=0:
					bit.append(str(int(snp.S&(2**x)!=0)) )
				else:
					bit.append('')
			ret.append(','.join(map(str, [scaffold, snp.bp, name]+[bR1]+[bR2]+[len(self.list)]+[snp.A, snp.B]+bit+['\n']) ) )
		return ''.join(ret)
	def get_state(self):
		bit=[]
		for x in range(0, self.numlines):
			if (~(self.M)&(2**x))!=0:
				if(int(self.S&(2**x)!=0)):
					bit.append('A')
				else:
					bit.append('B')
			else:
				bit.append('U')
		return ('\t'.join(bit) )


	def first(self, line):
		L=2**line
		for snp in self.list:
			if (~(snp.M)&L!=0):
				return snp.bp              
	def last(self, line):
		L=2**line
		for snp in reversed(self.list):
			if (~(snp.M)&L!=0):
				return snp.bp
	def get_index(self, line):
		return self.index[line]
	def set_index(self,line,bp):
		self.index[line]=bp
	def name(self):
		return self.scaffold+"_"+str(self.list[0].bp)+"_"+str(self.list[-1].bp)
 
		

def write_end(scaffold, group, numlines):
	ret=[]
	bit=[]
	AllBad=False
	for event in group:
		if event.bad:
			AllBad=True
#	if len(group)>5:
#		return ''
#	if AllBad:
#		return ''
	for event in group:
	#	event=event[0]
		bR1=str(event.start)
		bR2=str(event.stop)
		bp=event.stop
		if (len(group)%2==0):
			    TYPE="GC"
		else:
			    TYPE="CO"
		if event.bad:
			name=TYPE+"("+str(len(group) )+")"+str(event.line)+":#"+str(event.event)+"_BAD-"+str(event.count)
		else:
			name=TYPE+"("+str(len(group) )+")"+str(event.line)+":#"+str(event.event)+"_GOOD-"+str(event.count)
		bit=[]
		for x in range(0, numlines):
			if x==event.line:
				bit.append("1")
			else:
				bit.append('')
		ret.append(','.join(map(str, [scaffold, bp, TYPE]+[bR1]+[bR2]+[str(event.fivepSNPs)+"/"+str(event.threepSNPs)]+[name, 'Rec']+bit) ) )
	if ret!=[]:
		return '\n'.join(ret)+'\n'
	else: 
		return ''


def join(SNP1, SNP2):
	ret=SNP_BLOCK(SNP1.setbits, SNP1.numlines, SNP1.scaffold)
	[cR1, cR2, R1, R2]=getR(SNP1, SNP2)
	if cR1==0:
		ret.M=SNP1.setbits&( ~( (~SNP1.M)|(~SNP2.M) ) )
		ret.S=SNP1.setbits&( (SNP1.S & SNP2.M) | (~SNP2.M & SNP2.S) )
		ret.list=SNP2.list+SNP1.list
		return ret
	if cR2==0:
		SNP2.flip()
		ret.M=SNP1.setbits&( ~( (~SNP1.M)|(~SNP2.M) ) )
		ret.S=SNP1.setbits&( (SNP1.S & SNP2.M) | (~SNP2.M & SNP2.S) )
		ret.list=SNP2.list+SNP1.list
		return ret
		
numlines=len(SpermFile[0].split()[3:])
header=[]

for x in range(0, numlines):
	header.append("sp"+str(x+1) )

File=open("outR.csv", "w")
File.write("name,bp,type,# upstream recombinations,# downstream recombinations,#snp in group, allele1, allele2,"+','.join(header)+",\n")
File.close()

setbits=(2**numlines)-1

this=""

X=0

loci=[]

counts=[0]*numlines

while len(SpermFile)>0:
	#Read a line from the file.
	line=SpermFile.pop(0)
	line=line.split()
	
	if len(line)>1:
	#If the line is longer than 1, do something.
		try:
			calls=line[3:]
			for x in range(0, numlines):
				if (calls[x]!='-'):
					counts[x]+=1
			scaffold=line[0]
			bp=int(line[1])
			if scaffold==this:
				snp=SNP(bp, calls, setbits)
				if( (snp_block.add_snp(snp)) ):
					blocks.append(snp_block)
					snp_block=SNP_BLOCK(setbits, numlines, this)
					snp_block.add_snp(snp)
		except:
			print "Error!"
			print line

	if len(SpermFile)==0 or scaffold!=this:
		if this!="":
			blocks.append(snp_block)
			File=open("outR.csv", "a")
			if N50!=0:
				try:
					run=(int(this.split('_')[1] )<N50)
				except:
					print run
					print this.split('_')
					print "try running with the N50 0 option."
					quit()
			else:
				run=True
			
			if (run):
				#our first step is to recognize regions where recombination break points are roughly
				#similar on both the right and left region of the recombination block . . .

				bs=blocks[:]
				masked=[]
				bs.sort(key=lambda tup: len(tup.list), reverse=True)
				B=bs[0]

				#cut out poorly phased regions and bad sperm.

				#TODO make this user set
					
				not_done=True
				ML=0
				#Mask lines with low coverage
				avg=float(sum(counts))/float(numlines)
				if not (keep):
					for x in range(0, numlines):
						if (counts[x]<avg/4.0):
							print "masking ", x, " on scaffold ", this
						ML+=(counts[x]<avg/4.0)*2**x 
				setbits=(2**numlines)-1
				setbits=setbits&(~ML)

				while not_done:
					new_blocks=[]
					not_done=False
					for x in range(0, len(blocks) ):
						T=blocks[x]
						try:
							N=blocks[x+1]
						except:
							N=blocks[x]
						try:
							L=blocks[x-1]
						except:
							L=blocks[x]
						[cL1, cL2, L1, L2]=getR(L, T)
						[cN1, cN2, N1, N2]=getR(N, T)

						if min(cL1,cL2)==0 and (x>0) and len(new_blocks)>0:
							new_blocks.pop(-1)
							new_blocks.append( join(blocks[x], blocks[x-1]) )
							not_done=True
						else:
							if ( (len(T.list)<snp_thresh) and ( ( min(cL1, cL2)>max_rec) or (min(cN1, cN2)>max_rec)  ) ):
								print "appending"
								T.masked=True
								masked.append(T)
								not_done=True
							else:
								new_blocks.append(T)
					blocks=new_blocks

				print len(masked)
				masked.sort(key=lambda tup: tup.list[0].bp)

				#set all blocks in blocks to the state closest to (B)iggest.					
				for x in range(0, len(blocks) ):
					T=blocks[x]
					[cR1, cR2, R1, R2]=getR(B,T)
					if (cR2<cR1):
						blocks[x]=blocks[x].flip()
				
				#set all blocks in masked to the state closest to (B)iggest.
				for x in range(0, len(masked) ):
					T=masked[x]
					[cR1, cR2, R1, R2]=getR(B,T)
					if (cR2<cR1):
						masked[x]=masked[x].flip()

				#If everything is masked then there is nothing more to do.
				if len(blocks)==0:
					continue

				L=blocks[0]
				w=0
				A=SNP_BLOCK(setbits, numlines, this)

				hapblock=[]	
				for x in range(0, numlines):
					hapblock.append([])	
					hapblock[x].append(recombination(x, 0, 0, 0, 96) )
				#Phase locally. We are going to go through the blocks and put them in the state closest to their neighbor
				for x in range(0, len(blocks) ):
					T=blocks[x]
					try:
						N=blocks[x+1]
					except:
						N=blocks[x]

					[cR1, cR2, R1, R2]=getR(L, T)
					
					if len(blocks[x].list)<snp_thresh:
						#If the block has fewer markers in it than snp_thresh, don't use it to phasing down stream blocks,
						#but we still want to try to put it in the right phase...
						if (cR2<cR1):
							blocks[x]=T.flip()
							C1=R2
							cR1=cR2
					else:
						#otherwise phase down stream blocks based on this marker (i.e. set (L)ast to (T)his marker).
						L=T

					[cA1, cA2, A1, A2]=getR(A, T)
					A.M=A.setbits&( ~( (~A.M)|(~T.M) ) )
					A.S=A.setbits&( (A.S&T.M) | (~T.M&T.S) )

					#check lines for recombination.
					for y in range (0, numlines):
						if (2**y)&setbits !=0:
						#if the line isn't masked.
							if( (A1&2**y)!=0):
								hapblock[y].append(recombination(y, hapblock[y][-1].end, T.first(y), hapblock[y][-1].threepSNPs, max(min(cR1, cR2), min(nR1, nR2) ) ) )
							elif ((~T.M)&2**y)!=0:
								hapblock[y][-1].end=T.last(y)
								hapblock[y][-1].threepSNPs+=T.count[y]

					[nR1, nR2, R1, R2]=getR(N, L)
					
					File.write(T.write(this, "HAP", min(cR1, cR2), min(nR1, nR2) ) )

					if(bitCount(T.setbits&(~T.M))>minphased):
						loci.append(T)

					if (w<len(masked) ):
						print "Writing Masked", masked[w].list[0], N.list[0].bp
						while masked[w].list[0].bp<N.list[0].bp:
							bad=masked[w]
							[cR1, cR2, R1, R2]=getR(L, bad)
							[nR1, nR2, R1, R2]=getR(N, bad)
							File.write(masked[w].write(this, "BAD", min(cR1, cR2), min(nR1, nR2) ) )
							w+=1
							if w==len(masked):
								break
				for x in range(0, numlines):
					for y in range(x+1, numlines):
						for X in range(0, len(hapblock[x]) ):
							for Y in range(0, len(hapblock[y]) ):
								if (overlap(hapblock[x][X], hapblock[y][Y] ) ):
									hapblock[x][X].count+=1
									hapblock[y][Y].count+=1
				good=[]
				cluster=[]
				for x in range(0, numlines):
					good.append([])
					cluster.append([])
					
				for x in range(0, numlines):
					for event in hapblock[x]:
						if event.start==0:
							event.bad=True
							continue
						if error(event, 0.0000001, numlines, 0.01):
							event.bad=True
							good[event.line].append(event)
						else:
							good[event.line].append(event)

				for x in range(0, numlines):
					newcluster=[]
					laststop=0
					count=0
					for event in good[x]:
			#			if not (event.bad):
							if event.start-laststop<CGSIZE:
			#					if (event.?)?
								newcluster.append(event)
								laststop=event.stop
								event.event=count
							else:
								count+=1
			#					if (event.?)?
								cluster[x].append(newcluster[:])
								newcluster=[]
								newcluster.append(event)
								laststop=event.stop
								event.event=count
					cluster[x].append(newcluster[:])
				
				for x in range(0, numlines):
					for group in cluster[x]:
						File.write(write_end(this, group, numlines) )				
				File.close()
		counts=[0]*numlines
		snp=SNP(bp, calls, setbits)
		this=scaffold
		blocks=[]
		snp_block=SNP_BLOCK(setbits, numlines, this)
		snp_block.add_snp(snp)

LG={}
lg=0
LGlist={}
to_alloc=[]
MATERNAL={}
locidic={}
for locus in loci:
	locidic[locus.name()]=locus

for x in range(0, len(loci) ):
	for y in range(x+1, len(loci) ):
		L1=loci[x]
		L2=loci[y]
		[cR1, cR2, R1, R2]=getR(L1, L2)
		if min(cR1, cR2)<group_thresh and max(cR1,cR2)>numlines-group_thresh:
			if L1.name() in LG.keys():
				if L2.name() in LG.keys():
					from_lg=LG[L1.name()]
					to_lg=LG[L2.name()]
					if from_lg==to_lg:
						continue
					for key in LGlist[from_lg]:
						LG[key]=to_lg
						LGlist[to_lg].append(key)
					if(cR2<cR1):
						for key in LGlist[from_lg]:
							locidic[key].flip()
					del LGlist[from_lg]
					to_alloc.append(from_lg)
				else:
					LG[L2.name()]=LG[L1.name()]
					LGlist[LG[L2.name()]].append(L2.name())
					if(cR2<cR1):
						L2.flip()
			elif L2.name() in LG.keys():
				LG[L1.name()]=LG[L2.name()]
				LGlist[LG[L2.name()]].append(L1.name())
				if(cR2<cR1):
					L1.flip()
			else:
				try:
					to_lg=to_alloc.pop()
				except:
					to_lg=lg
					lg+=1
				LGlist[to_lg]=[L2.name(), L1.name()]
				LG[L2.name()]=to_lg
				LG[L1.name()]=LG[L2.name()]
				if(cR2<cR1):
					L2.flip()
print "found "+str(len(LGlist))+" tenetive linkage groups."
File=open("MST-data.txt", 'w')
File.write("population_type DH\npopulation_name LG\ndistance_function kosambi\ncut_off_p_value 0.000001\nno_map_dist 15.0\nno_map_size 0\nmissing_threshold 1.00\nestimation_before_clustering no\ndetect_bad_data yes\nobjective_function COUNT\nnumber_of_loci "+str(len(loci))+"\nnumber_of_individual "+str(numlines)+"\n")
File.write("\nlocus_name\t"+'\t'.join(header)+'\n')

for x in range(0, len(loci) ):
	File.write(loci[x].name()+"\t"+loci[x].get_state()+'\n')


File.close()
