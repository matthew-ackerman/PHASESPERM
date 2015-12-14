
pdf("outAll.pdf")

data<-read.csv("outR.csv")

SI<-6
scaffold<-data[,1]
pos<-data[,2]/1000
type<-data[,3]
calls<-data[,9:ncol(data)]

A<-0.5
L<-ncol(data)-9
N<-10
OFF<-10

x<-min(pos)
y<-max(pos)


#iterate over the scaffolds
for (s in unique(scaffold) ){

#
for (j in 0:10){
	plot(pos[scaffold==s & type=="EMPTY"],calls[,1][scaffold==s & type=="EMPTY"], pch='|', bty='n', xlim=c(x, y), ylim=c(-(j+1)*N*SI+0.25*OFF, -j*N*SI-4.75+0.25*OFF ), xlab="position (kb)", ylab="sperm", yaxt='n', main=s)


#for (z in 1:L){
#	rect(pos[scaffold==s & type=="ERR"], pos[scaffold==s & type=="ERR"], ,col="red")
#}

for (z in 1:L){
	rect(0,0.5-z*SI, max(pos[scaffold==s]), -1.5-z*SI, col="gray")
}

for (z in 1:L){
	points(pos[scaffold==s & type=="BAD"], -calls[,z][scaffold==s & type=="BAD"]-z*SI, col=ifelse( (calls[,z][scaffold==s & type=="BAD"])==0, "gray58", "gray58") , cex=A, pch='|')
	points(pos[scaffold==s & type=="HAP"], -calls[,z][scaffold==s & type=="HAP"]-z*SI, col=ifelse( (calls[,z][scaffold==s & type=="HAP"])==0, "black", "red") , cex=A, pch='|')
	points(pos[scaffold==s & type=="REC"], -calls[,z][scaffold==s & type=="REC"]-z*SI+2, col="green", cex=A, pch=6)
	points(pos[scaffold==s & type=="GEC"], -calls[,z][scaffold==s & type=="GEC"]-z*SI+2, col="yellow", cex=A, pch=6)
	points(pos[scaffold==s & type=="ERR"], -calls[,z][scaffold==s & type=="ERR"]-z*SI+2, col="red", cex=A, pch=6)
}

for (z in 0:L){
	text( x/1000-(x-y)*0.01, -0.45-SI*z, z, cex=0.75)

}
}
}
warnings()
dev.off()
