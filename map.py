import math
import sys

File=open(sys.argv[1])
calls=[]

def bitCount(int_type):
	count = 0
	while(int_type):
		int_type &= int_type - 1
		count += 1
	return(count)

class call:
	def __init__ (self, scaffold, start, stop, ID, MASK, SET):
		self.scaffold=scaffold
		self.start=int(start)
		self.stop=int(stop)
		self.ID=int(ID)
		self.MASK=int(MASK)
		self.SET=int(SET)
	
def get_LOD(A, B):
	R=min(bitCount(A.SET&( ( (~A.MASK)&(~B.MASK) )&(B.ID^A.ID) ) ), bitCount(A.SET&( ((~A.MASK)&(~B.MASK))&(B.ID^(~A.ID) ) ) ) )
	C=bitCount( A.SET&(~A.MASK)&(~B.MASK) )
	Theta=float(R)/float(C)
	if Theta==0:
		return [-0.5*math.log(1-2*Theta)*100, 1-math.log(0.25, 10)*C]
	elif Theta==0.5:
		return ["NaN", math.log(Theta/2.0, 10)*R+math.log((1-Theta)/2.0, 10)*(C-R)-math.log(0.25, 10)*C]
	else:
		return [-0.5*math.log(1-2*Theta)*100, math.log(Theta/2.0, 10)*R+math.log((1-Theta)/2.0, 10)*(C-R)-math.log(0.25, 10)*C]

for line in File:
	line=line.split(',')
	try:
		calls.append(call(line[0], line[1], line[2], line[3], line[4], line[5]) )
	except:
		print line
	
for x in range(0, len(calls) ):
	A=calls[x]
	for y in range(x+1, len(calls) ):
		B=calls[y]
		LOD=get_LOD(A, B)
		if LOD[1]>4:
			print ','.join(map(str, [A.scaffold, A.start, B.scaffold, B.start]+LOD) )
