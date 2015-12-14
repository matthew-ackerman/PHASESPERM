import sys

File=open(sys.argv[1])

out=[]

lastlg=0
for line in File:
	line=line.strip('\n').split()
	if len(line)>1:
		if line[0]=="group":
			LG=line[1]
		elif line[0][0:2]=="sc":
			scf=line[0].split('_')[0]
			start=int(line[0].split('_')[1])
			stop=int(line[0].split('_')[2])
			cm=float(line[1])
			bpcode=start+10**8*int(scf[2:])
			lgcode=cm+10**5*int(LG[2:])
			if lastlg!=lgcode:
				out.append([LG, scf, start, stop, cm, 0, bpcode, lgcode] )
			lastlg=lgcode
out=sorted(out,key=lambda x: x[6])

for x in range(1, len(out) ):
	out[x][5]=out[x-1][5]+1

out=sorted(out,key=lambda x: x[7])

start_bp=out[0][2]
stop_bp=out[0][3]
start_cm=out[0][4]
forward=False
reverse=False
for x in range(1, len(out) ):
	#print out[x]
	if (forward):
		if out[x][5]!=(out[x-1][5]+1):
			print ','.join(map(str, [out[x][0], out[x][1], start_bp, out[x][3], start_cm, out[x][4]]))
			start_bp=out[x][2]
			stop_bp=out[x][3]
			start_cm=out[x][4]
			foward=False
	elif (reverse): 
		if out[x][5]!=(out[x-1][5]-1):
			print ','.join(map(str, [out[x][0], out[x][1], start_bp, out[x][2], start_cm, out[x][4]]))
			start_bp=out[x][2]
			stop_bp=out[x][3]
			start_cm=out[x][4]
			foward=False
	else:
		if out[x][5]==(out[x-1][5]-1):
			reverse=True
		if out[x][5]==(out[x-1][5]-1):
			forward=True
		else:
			print ','.join(map(str, [out[x][0], out[x][1], start_bp, stop_bp, start_cm, start_cm]))
			start_bp=out[x][2]
			stop_bp=out[x][3]
			start_cm=out[x][4]
