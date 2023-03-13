rm(list=ls())
library(tidyr)
library(dplyr)

setwd("C:/Users/USUARIO/Desktop/FPIES/FPIES-IBS")

# Lectura de archivos

anot <- read.table("/path/to/file/anotation.txt",header = TRUE, fill = TRUE)
bim <- read.table("file.bim", header = FALSE, 
                  col.names = c("CHR", "SNP", "CT", "POS", "REF", "ALT"))

# A la hora de unir los dos ficheros, se genera un código chr-pos para evitar errores

code_bim <- bim[,c("CHR", "POS")]
code_bim <-unite(code_bim, code_bim,c(1:2), sep="-", remove=TRUE)
bim <- cbind(bim, code_bim)

code_anot <- anot[,c("CHR", "POS")]
code_anot[1] <- gsub("chr", "", anot$CHR)
code_anot <-unite(code_anot, code_anot,c(1:2), sep="-", remove=TRUE)
anot <- cbind(anot, code_anot)

# Unión de los dos ficheros a través de ese código

file_merge <- merge(bim, anot, by.x = "code_bim", by.y = "code_anot")

# Eliminar variantes duplicadas, indels, trialélicos, etc

file_merge_acgt <- file_merge[(file_merge$ALT.x%in%c("A", "T", "C", "G") & 
                                 (file_merge$REF.x%in%c("A", "T", "C", "G"))),]
sum(duplicated(file_merge$SNP))
dupe <- file_merge_acgt[duplicated(file_merge_acgt$SNP),]
file_merge_nodup <- file_merge_acgt[!(file_merge_acgt$SNP%in%dupe$SNP),]

data_filtered <- file_merge_nodup[,c("CHR.x", "SNP", "POS.x", "ID", "GeneName")]

# Dividir la columna GeneName por "|" para obtener todos los genes

genenames <- strsplit(data_filtered$GeneName, split = "|", fixed = TRUE)

# Se utiliza la función lapply que selecciona cada elemento de la lista y le aplica
# la función length

long <- unlist(lapply(genenames, length))
genenames_str <- unlist(genenames)

# Lista de SNPs

snps <- data_filtered$SNP

# Se repite el snp tantas como veces como genes haya para esa variante

snps_list <- cbind(rep(snps, long))

snps_genes <- cbind(snps_list, genenames_str)
snps_genes <- data.frame(snps_genes)

# Filtrar aquellas filas que contienen NAs

data_depured <- snps_genes[which(snps_genes$genenames_str!="NA"),]
data_depured <- data_depured[which(data_depured$genenames_str!="."),]

# Generar código snp-gen 

data_depured$code <- paste(data_depured$V1, data_depured$genenames_str, sep = "/")

# Eliminar duplicados (en algunas variantes se puede repetir el gen en el que aparece)

sum(duplicated(data_depured$code))
data_nodupe <- data_depured[!duplicated(data_depured$code),]

# Eliminar columna de código

data_nodupe <- data_nodupe[,c("V1", "genenames_str")]
data_nodupe[,3] <- data_nodupe[,1]
data_nodupe[,1] <- data_nodupe[,2]
data_nodupe[,2] <- data_nodupe[,3]

# Ordenar alfabéticamente los genes

setid <- data_nodupe[order(data_nodupe$V1),]

# Se debe guardar primero los genes y luego los SNPs

setid<- setid[,c("V1", "genenames_str")]

# Escribir fichero SetID

write.table(setid, "file.SetID", col.names = FALSE, row.names = FALSE, quote = FALSE)

