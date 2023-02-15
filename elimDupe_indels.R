rm(list=ls())

library(dplyr)
library(tidyr)
library(stringr)

setwd("C:/Users/usuario/Desktop/IDIS/FPIES")

# Lectura de archivos

bim2 <- read.table("kfe_fpies_depured.bim", header=FALSE, col.names=c("CHR", "VAR", "CENT", "POS", "ALT", "REF"))

# Crear código para cada variante de la forma CHR-POS

code_bim <- bim[,c("CHR", "POS")]
code_bim <-unite(code_bim, code_bim,c(1:2), sep="-", remove=TRUE)
bim <- cbind(bim, code_bim)

# Encontrar duplicados en la base de datos

dupe <- code_bim[duplicated(code_bim$code),]
dupe <- data.frame(dupe)
colnames(dupe) <- "code"
colnames(code_fpies) <- "code"

#  Detectar variantes con alelos distintos a ACGT

rs <- bim[(bim$code%in%dupe$code),]
indels <- bim[!(bim$ALT%in%c("A", "T", "C", "G")& bim$REF%in%c("A", "T", "C", "G")),]
indels <- bind_rows(indels, rs)
indels <- indels[-7]

# Escribir fichero con las variantes que se quieren eliminar

write.table(chr_rm, file="snps_dupe_indels.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)

