rm(list = ls())
library(scales)

setwd("C:/Users/USUARIO/Desktop")

## Lectura de archivos, aqui no tengo poblaciones

samples <-read.csv("FPIES/KFE/khor_sample_sheet_2.csv",sep = ";",header=TRUE)
fpies <- read.table("FPIES/KFE/Exomas_BioFPIES.txt",header = TRUE, fill = TRUE)
ibs<-read.table("FPIES/HC_FPIES/mds_hck.mibs")
id<-read.table("FPIES/HC_FPIES/mds_hck.mibs.id");colnames(id)<-c("FID","IID")


#inf <- inf[match(id[,2], inf[,1]),]
#inf[which(is.na(inf[,2])),]<-cbind(id[which(is.na(inf[,2])),1],"CON")

samples <- samples[,c("Khor", "Code", "Cohort")]
fpies <- fpies[,c("BIOBANK_ID", "EXOME", "COHORT")]
colnames(samples) <- c("ID", "Code", "Cohort")
colnames(fpies) <- c("ID", "Code", "Cohort")

fpies[1] <- fpies[2] 

samples_conj <- rbind(fpies, samples)

inf<-samples_conj[match(id[,"FID"],samples_conj[,1]),]
inf_fpies <-match(fpies[,2],id[,"FID"], )
inf_khor <- match( samples[,1],id[,"FID"],)

inf <-inf[which(inf$ID!="NA"),]

##########################
Cobertura<-sample(c("NA"), size = 42, replace = TRUE)
inf <- cbind(inf, Cobertura)
low_cobert <- inf[which(startsWith(inf$Code, "B") | startsWith(inf$Code, "A")),]
inf$Cobertura <- ifelse(inf$Code%in%low_cobert$Code, "Baja cobertura - enterocolitis", "Alta cobertura - enterocolitis")
#########################

inf$Cohort<-ifelse(inf$Code%in%fpies$ID, "FPIES-A", "Control")
## Vector con los tipos de cohortes
ls<-unique(inf[,3])

pch<-rep(1:8,length.out=length(ls))
pch<-1:length(ls);col<-rainbow(length(ls))


## Se calcula la matriz de distancias
m<-match(inf[,3],ls)
mds<-cmdscale(as.dist(1-ibs),k=nrow(ibs)-1,add=T,eig=T)
var_ex <-100*mds$eig/sum(mds$eig)

pdf("FPIES/HC_FPIES/mds_prueba1.pdf",width=15,height=7*20/18)
layout(rbind(c(1,1),c(2,3)),height=c(1,9))
par(mar=rep(0,4));plot.new()
legend("left",bty="n",col=col,pch=pch,legend=ls,nc=4,cex=1)
par(mar=c(5,5,1,1))
plot(mds$points[,1:2],col=alpha(col[m],0.35),pch=pch[m],xlab=paste("DIMENSION 1"," (",round(var_ex[1],2),"%)",sep=""),ylab=paste("DIMENSION 2"," (",round(var_ex[2],2),"%)",sep=""),cex.lab=1.2,cex.axis=1.2,lwd=1.5,cex=1.5)
dev.off()

