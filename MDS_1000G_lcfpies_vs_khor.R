rm(list = ls())
library(scales)

setwd("C:/Users/USUARIO/Desktop/FPIES")

#samples <-read.csv("Final_Samples_Singap_khor.csv",sep = ";",header=TRUE)
#codes <-read.csv("khor_codes.csv", sep = ";", header=TRUE)

samples <-read.csv("KFE/khor_sample_sheet_2.csv",sep = ";",header=TRUE)
fpies <- read.table("KFE/Exomas_BioFPIES.txt",header = TRUE, fill = TRUE)
info <- read.table("LC_FPIES/lcfpies_1000g/Informacion_1000G.txt", header =TRUE, fill=TRUE, sep = "\t", 
                  colClass = "character")

ibs<-read.table("Total_FPIES/mds_fpies_khor_1000g.mibs")
id<-read.table("Total_FPIES/mds_fpies_khor_1000g.mibs.id");colnames(id)<-c("FID","IID")

samples <- samples[,c("Khor", "Code", "Cohort")]
fpies <- fpies[,c("BIOBANK_ID", "EXOME", "COHORT")]

colnames(samples) <- c("ID", "Code", "Cohort")
colnames(fpies) <- c("ID", "Code", "Cohort")
low_cobert <- fpies[which(startsWith(fpies$Code, "18") | startsWith(fpies$Code, "19")),]

samples_conj <- rbind(low_cobert, samples)


colnames(info)[2] <- "IID"

inf <- info[,c("IID", "Population")]
inf <- inf[match(id[,2], inf[,1]),]

## Se diferencian los Controles (KHOR) y casos (FPIES)

#inf[which(is.na(inf[,2])),]<-cbind(id[which(is.na(inf[,2])),1],"CON")

inf[which(is.na(inf[,2])),]<-cbind(samples_conj[which(is.na(inf[,2])),1])
inf[,2] <- ifelse(inf[,2]%in%low_cobert$ID, "CASOS", inf[,2])
inf[,2] <- ifelse(inf[,2]%in%samples[,1], "CON", inf[,2])


ls<-unique(inf[,2])

#Identifico la etiqueta casos
w <- which(ls == "CASOS")

pch<-rep(1:18,length.out=length(ls));col<-rainbow(length(ls))

#Asigno a casos un color y simbolo concreto
pch[w]=16; col[w]=1

## Se calcula la matriz de distancias
m<-match(inf[,2],ls)
mds<-cmdscale(as.dist(1-ibs),k=nrow(ibs)-1,add=T,eig=T)
var_ex <-100*mds$eig/sum(mds$eig)

library("scales")

setwd("C:/Users/usuario/Desktop/FPIES")
pdf("HC_FPIES/hcfpies_1000g/1000g_hcfpies_khor_selectpobs.pdf",width=15,height=7*20/18)
layout(rbind(c(1,1),c(2,3)),height=c(1,9))
par(mar=rep(0,4));plot.new()
legend("center",bty="n",col=col,pch=pch,legend=ls,nc=12,cex=1)
par(mar=c(5,5,1,1))
w<-which(inf[,2]=="CASOS")
plot(mds$points[,1:2],col=alpha(col[m],0.35),pch=pch[m],xlab=paste("DIMENSION 1"," (",round(var_ex[1],2),"%)",sep=""),ylab=paste("DIMENSION 2"," (",round(var_ex[2],2),"%)",sep=""),cex.lab=1.2,cex.axis=1.2,lwd=1.5,cex=1.5)
points(mds$points[w,1:2],col=alpha(col[m[w]],0.35),pch=pch[m[w]],cex=1.5)

plot(mds$points[,3:4],col=alpha(col[m],0.35),pch=pch[m],xlab=paste("DIMENSION 3 (",round(var_ex[3],2),"%)",sep=""),ylab=paste("DIMENSION 4"," (",round(var_ex[4],2),"%)",sep=""),cex.lab=1.2,cex.axis=1.2,lwd=1.5,cex=1.5)
points(mds$points[w,3:4],col=alpha(col[m[w]],0.35),pch=pch[m[w]],cex=1.5)
dev.off()

#Con points despues de cada plot vuelvo a representar los puntos de casos para que queden por encima de
#la gráfica y se puedan apreciar mejor

