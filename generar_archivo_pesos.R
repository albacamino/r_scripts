rm(list=ls())
library(qqman)
library(dplyr)
library(tidyr)
library(stringr)
setwd("/file/to/directory")

# Lectura de archivos

assoc <- read.table("file.SetID", header = FALSE)
bim <- read.table("filet.bim", header = FALSE)
anot <- read.table("/path/to/file/anotation.txt", 
                   header = TRUE)

# Se crea un código para emparejar los dos archivosç

code_todos <- bim[,c("V1", "V4")]
code_todos <- unite(code_todos, code_todos, c(1:2), sep="-", remove = TRUE)
bim<- cbind(bim, code_todos)

code_anot <- anot[,c("CHR", "POS")]
code_anot[1] <- gsub("chr", "", anot$CHR)
code_anot <- unite(code_anot, code_anot, c(1:2), sep="-", remove = TRUE)
anot<- cbind(anot, code_anot)
colnames(anot)[12] <- "code"

# Se matchean los dos archivos por la columna código

common <- merge(anot,bim, by.x = "code", by.y = "code_todos" )
common_genes <- common[,c("V2", "AnnoType", "Consequence", "GeneName", "PHRED")]
colnames(common_genes)[1] <- "SNP"
colnames(assoc)[2] <- "SNP"

# Se comprueba con el fichero SetID

union<- merge(assoc, common_genes, by.x = "SNP", by.y = "SNP")

# Nos quedamos con las columnas de interés

weights <- union[,c("SNP", "PHRED")]

# Se escribe el fichero

write.table(weights, "file.txt", col.names = FALSE, 
            row.names = FALSE, sep = " ",quote = FALSE)

