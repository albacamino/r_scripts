rm(list=ls())
setwd("C:/Users/USUARIO/Desktop/FPIES")

# Lectura de archivos
inf<-read.table("LC_FPIES/lcfpies_1000g/Informacion_1000G.txt",sep="\t",head=T,colClass="character")[,c(2,7)]
fam<-read.table("Total_FPIES/fpies_khor_1000g_selectpobs.fam")
fam_klc <- read.table("Total_FPIES/fpies_vs_khor_nodup.fam",header = FALSE, 
                      col.names = c("FID","IID", "V3", "V4", "V5", "V6"))
samples_khor <-read.csv("KFE/khor_sample_sheet_2.csv",sep = ";",header=TRUE)
fpies <- fam_klc[which(startsWith(fam_klc[,1], "A") | startsWith(fam_klc[,1], "B") | startsWith(fam_klc[,1], "18") | startsWith(fam_klc[,1], "19")),]

# Match  entre el archivo .fam y el de informacion de 1000G
inf<-inf[match(fam[,2],inf[,1]),] 

# Se rellenan las filas con NA con los datos de casos y controles, que no aparecen
# en el archivo de información de 1000G
inf[which(is.na(inf[,2])),]<-cbind(fam_klc[which(is.na(inf[,2])),1])
inf[,2] <- ifelse(inf[,2]%in%fpies[,1], "CASOS", inf[,2])
inf[,2] <- ifelse(inf[,2]%in%samples_khor[,1], "CON", inf[,2])

inf_fpies <- inf[which(startsWith(inf[,1], "A") | startsWith(inf[,1], "B")),]

inf[which(fam[,ncol(fam)]==2),2]<-"CAP"
k<-4

# Se abre el archivo de ancestralidad de ADMIXTURE (el .Q tiene el mismo orden que el .fam)
q<-read.table("Total_FPIES/fpies_khor_1000g_selectpobs.4.Q")

ls<-unique(as.character(inf[,2]))
ls<-c("GBR", "IBS", "CEU", "CHS", "CDX", "JPT", "CLM", "PEL", "MXL", "ESN", "MSL", "YRI","CON","CASOS" ) 
inf<-inf[order(match(inf[,2],ls)),]

# Se ordenan los datos en función del orden de las pobs y se le pasa al .fam y al .Q para mantener el orden
order<-match(inf[,1],fam[,1])
order_fpies <- match(inf[,1], fpies[,1])
order_lc<- order_fpies[!is.na(order_fpies)]

fam<-fam[order,]
col<-rainbow(4)
q<-q[order,]



# Matriz de gráficos en la figura
mat <- matrix(c(1,1,2), ncol = 3)

layout(mat = mat)

# PRIMER GRÁFICO
# Promedio de cada de una de las ancestralidades en cada una de las componentes
frequencies<-array(dim=c(length(ls),k))
rownames(frequencies)=ls
for(i in 1:length(ls)){
  a<-q[which(inf[,2]==ls[i]),]
  frequencies[i,]<-colMeans(a)
}
freq<-order(apply(frequencies,2,which.max));q<-q[,freq]
h1<-as.numeric()

# Ordenar indvs dentro de cada pob e introducir separadores entre las poblaciones
for(i in 1:(length(ls)-1)){
  a<-q[which(inf[,2]==ls[i]),]
  h1<-rbind(h1,cbind(a[rev(do.call("order",a)),],rep(0,nrow(a))),c(rep(0,ncol(a)),1),
            c(rep(0,ncol(a)),1))
}
a<-q[which(inf[,2]==ls[length(ls)]),]
# Los casos son los ultimos en el dataframe asi que se guardan en b
b<-cbind(a[rev(do.call("order",a)),],rep(0,nrow(a)))


par(mar=c(4,5,1.5,0.00005))
mp<-barplot(100*t(h1),col=c(col,1),space=0,border=NA,xaxt='n',ylab="ANCESTRY (%)",
            cex.axis=1.25,cex.lab=1.25)
a<-1:length(ls)-1

for(i in 1:length(a)){
  a[i]<-mean(mp[which(inf[,2]==ls[i])+2*(i-1)])
}
text(a, par("usr")[3] - 4, srt = 90, adj = 1, labels = ls, xpd = TRUE,cex=0.85)

par(mar=c(4,0.0005,1.5,1.5))


mp<-barplot(100*t(b),col=c(col,1),space=0,border=NA,xaxt='n',cex.axis=1.25,cex.lab=1.25,
            axes = FALSE, las = 2)
text(x = 10 , y = par("usr")[3]- 4, labels = "CASOS", cex = 0.85, xpd = TRUE, srt = 90, adj = 1)


