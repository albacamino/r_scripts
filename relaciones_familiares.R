rm(list = ls())
library("Ternary")
library("scales")

setwd("C:/Users/USUARIO/Desktop/FPIES")
genome <- read.table("C:/Users/USUARIO/Desktop/FPIES/LC_FPIES/lcfpies_genome.genome", header=TRUE)

im <-rowSums(genome[,c("HOMHOM","HETHET")])
a <-genome[,"HETHET"]/im
v <-round(1-pnorm((a-2/3)/sqrt(2/(9*im))),3)
w0 <-which(v<0.05/length(v))

pdf("C:/Users/USUARIO/Desktop/FPIES/LC_FPIES/relaciones_fam_LC_fpies.pdf")
TernaryPlot(axis.labels = seq(0, 1, by = .1),alab = 'IBD0', blab = 'IBD1', clab = 'IBD2')
TernaryPoints(apply(genome[,c("Z0","Z1","Z2")],2,as.numeric),col=alpha("red",0.5),pch=16)

b<-cbind(c(0,0,1,0.25,0.5,0.75,0.5625,0.4375,0.65625),c(0,1,0,0.5,0.5,0.25,0.375,0.5,0.3125))
TernaryPoints(cbind(b,1-rowSums(b)),pch=c(15:18,1:5),cex=1.5)
TernaryPoints(genome[,c("Z0","Z1","Z2")],col=alpha("red",0.5),pch=16)
legend("topright",pch=c(15:18,1:5),bty="n",	legend=c("Duplicated","Parent-son","Unrelated","Siblings","Avuncular"," First cousins","Double first cousins",	"Half-siblings cousins","First-second cousins"),pt.cex=1.2,y.intersp=1.4,cex=0.75)
dev.off()

pdf("C:/Users/USUARIO/Desktop/FPIES/LC_FPIES/Kinship_plot_LC_fpies.pdf",width=7,height=7)
plot(0,0,type="n",xlab = expression(kappa[0]), ylab = expression(kappa[1]),main="",frame.plot=FALSE,	xlim=c(0,1),ylim=c(0,1),cex.lab=1.3)

segments(0,0,1,0,lwd=2);segments(0,0,0,1,lwd=2);segments(1,0,0,1,lwd=2)
points(genome[-w0,c("Z0","Z1")],col=alpha("red",0.5),pch=16,cex=1.5)

points(c(0,0,0.25,0.5,0.75,0.5625,0.4375,0.65625,1),c(0,1,0.5,0.5,0.25,0.375,0.5,0.3125,0),pch=c(15:18,1:5),cex=1.8)
legend("topright",pch=c(15:18,1:5),bty="n",legend=c("Duplicated","Parent-son","Siblings","Avuncular"," First cousins","Double first cousins",	"Half-siblings cousins","First-second cousins","Unrelated"),pt.cex=1.8,y.intersp=1.4)

k0<-seq(0,1,by=0.01)
k1.1<-2*(-k0+sqrt(k0))
k1.2<-2*(-k0-sqrt(k0))
lines(k0,k1.1,lwd=2,lty=2)
points(genome[w0,c("Z0","Z1")],col=alpha("red",0.5),pch=16,cex=1.5)
dev.off()

