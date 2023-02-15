rm(list=ls())
library(dplyr)
library(tidyr)
library(stringr)
setwd("C:/Users/USUARIO/Desktop/FPIES/KFE")

# Lectura de archivos

bim <- read.table("kfe-justsnp_3.bim", header = FALSE, col.names=c("CHR", "VAR", "CENT", "POS", "ALT", "REF"))

# Crear código para cada variante en el archivo

code_bim <- bim[,c("CHR", "POS")]
code_bim <-unite(code_bim, code_bim,c(1:2), sep="-", remove=TRUE)
bim <- cbind(bim, code_bim)

# Variantes duplicadas

dupe <- code_bim[duplicated(code_bim$code),]
dupe <- data.frame(dupe)
colnames(dupe) <- "code"
colnames(code_fpies) <- "code"

# Detectar variantes con alelos distintos a ACGT

rs <- bim[(bim$code%in%dupe$code),]
indels <- bim[!(bim$ALT%in%c("A", "T", "C", "G")& bim$REF%in%c("A", "T", "C", "G")),]
indels <- bind_rows(indels, rs)
indels <- indels[-7]

indels <- indels[2]

write.table(indels, "indels_remove.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)


