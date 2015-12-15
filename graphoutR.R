
pdf("outAll.pdf")

data<-read.csv("outR.csv")

SI<-6				#Scale Needs to be >4, since 2 for the '|' 1 for the triangles, and at least 1 for blank space.
scaffold<-data[,1]		#All the scaffold names
pos<-data[,2]/1000		#pos in kb.
type<-data[,3]			#Code for the type of data "HAP" for haplotype call (0/1) "BAD" for a bad site, "REC" for a CO
calls<-data[,9:ncol(data)]	#

A<-0.5				#Size of vertical marks
L<-ncol(data)-9			#Number of sperm
N<-10				#Sperm per page
OFF<-10				#?

x<-min(pos)			#?
y<-max(pos)			#?


#iterate over the scaffolds
for (s in unique(scaffold) ){
	for (j in 0:( (L-1)/N) ){
		print(s)
		print(j)
		plot(pos[scaffold==s & type=="EMPTY"],calls[,1][scaffold==s & type=="EMPTY"], pch='|', bty='n', xlim=c(x, y), ylim=c(-N*SI-0.25*OFF, 4.75+0.25*OFF ), xlab="position (kb)", ylab="sperm", yaxt='n', main=s)
#for (z in 1:L){
#	rect(pos[scaffold==s & type=="ERR"], pos[scaffold==s & type=="ERR"], ,col="red")
#}

		for (z in (j*N+1):min((j+1)*N, L) ){
			rect(0,0.5-(z-j*N)*SI, max(pos[scaffold==s]), -1.5-(z-j*N)*SI, col="gray")
		}
		for (z in (j*N+1):min((j+1)*N, L) ){
			print(z)
			points(pos[scaffold==s & type=="BAD"], -calls[,z][scaffold==s & type=="BAD"]-(z-j*N)*SI, col=ifelse( (calls[,z][scaffold==s & type=="BAD"])==0, "gray58", "gray58") , cex=A, pch='|')
			points(pos[scaffold==s & type=="HAP"], -calls[,z][scaffold==s & type=="HAP"]-(z-j*N)*SI, col=ifelse( (calls[,z][scaffold==s & type=="HAP"])==0, "black", "red") , cex=A, pch='|')
			points(pos[scaffold==s & type=="CO"], -calls[,z][scaffold==s & type=="CO"]-(z-j*N)*SI+2, col="green", cex=A, pch=6)
			points(pos[scaffold==s & type=="GC"], -calls[,z][scaffold==s & type=="GC"]-(z-j*N)*SI+2, col="orange", cex=A, pch=6)
		#	points(pos[scaffold==s & type=="ERR"], -calls[,z][scaffold==s & type=="ERR"]-z*SI+2, col="red", cex=A, pch=6)
		}

		for (z in (j*N+1):min((j+1)*N, L) ){
			text( x-(y-x)*0.025, -0.45-SI*(z-j*N), z, cex=0.75)
		}
	}
	break
}
warnings()
dev.off()
