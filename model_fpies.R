rm(list=ls())
library(qqman)
library(dplyr)
library(tidyr)
library(stringr)
setwd("C:/Users/USUARIO/Desktop/FPIES")

# Lectura de archivos

model <- read.table("Total_FPIES/HC_FPIES/hcfpies_assoc_hwe/model_hcfpies.model", header = TRUE)
bim <- read.table("Total_FPIES/HC_FPIES/hcfpies_assoc_hwe/hcfpies_vs_khor.bim", header = FALSE, 
                  col.names = c("CHR", "VAR", "CENT", "POS", "ALT", "REF"))
assoc <- read.table("Total_FPIES/HC_FPIES/hcfpies_assoc_hwe/hcfpies_khor_varselect_assocc.assoc", header = TRUE)
anot <- read.table("Total_FPIES/total_fpies_assoc/fpies_annot.txt", header = TRUE)


# Eliminar los NAs

model_depured <- model[complete.cases(model),]

# Añadir la posición de las variantes y el F_U para filtrar

model_bp <- merge(assoc, model_depured, by.x = "SNP", by.y = "SNP")
model_bp <- model_bp[,c("CHR.x", "SNP", "BP", "F_A", "F_U", "A1.x", "A2.x", "TEST", "AFF", "UNAFF", "CHISQ.x", "DF", "P.x" )]

colnames(model_bp)[1] <- "CHR"
colnames(model_bp)[6] <- "A1"
colnames(model_bp)[7] <- "A2"
colnames(model_bp)[11] <- "CHISQ"
colnames(model_bp)[13] <- "P"

# Separar en variables diferentes según el modelo aplicado

dom <- model_bp[which(model_bp$TEST=="DOM"),]
rec <- model_bp[which(model_bp$TEST=="REC"),]
geno <- model_bp[which(model_bp$TEST=="GENO"),]
allelic <- model_bp[which(model_bp$TEST=="ALLELIC"),]
trend <- model_bp[which(model_bp$TEST=="TREND"),]

# Filtrar por la columna F_U

dom_filter <- dom[which(dom$F_U>0.05),]
rec_filter <- rec[which(rec$F_U>0.05),]
geno_filter <- geno[which(geno$F_U>0.05),]
allelic_filter <- allelic[which(allelic$F_U>0.05),]
trend_filter <- trend[which(trend$F_U>0.05),]

# Manhattan plots con cada uno de los test 

manhattan(dom_filter, chr="CHR", bp="BP", snp="SNP", p="P", annotatePval = 0.01, 
          main = "Manhattan plot - FPIES ~ Dominant")
manhattan(rec_filter, chr="CHR", bp="BP", snp="SNP", p="P", annotatePval = 0.01, 
          main = "Manhattan plot - FPIES ~ Recessive")
manhattan(geno_filter, chr="CHR", bp="BP", snp="SNP", p="P", annotatePval = 0.01, 
          main = "Manhattan plot - FPIES ~ Genotypic")
manhattan(allelic_filter, chr="CHR", bp="BP", snp="SNP", p="P", annotatePval = 0.01, 
          main = "Manhattan plot - FPIES ~ Allelic")
manhattan(trend_filter, chr="CHR", bp="BP", snp="SNP", p="P", annotatePval = 0.01, 
          main = "Manhattan plot - FPIES ~ Cochran-Armitage trend")

# Seleccionar las variantes que superan la linea azul (suggestive line) en el test alélico

table_genes <- allelic_filter[which((-log10(allelic_filter$P))>(-log10(1e-05))),]

code_table <- table_genes[,c("CHR", "BP")]
code_table <- unite(code_table, code_table, c(1:2), sep="-", remove = TRUE)
table_genes <- cbind(table_genes, code_table)
colnames(table_genes)[14] <- "code"

# Con gsub busco el primer elemento de la función y lo sustituyo por el segundo elemento
# en el objeto, que se corresponde con el terce parámetro.

anot[1] <- gsub("chr", "", anot$CHR)
code_anot <- anot[,c("CHR", "POS")]
code_anot <- unite(code_anot, code_anot, c(1:2), sep="-", remove = TRUE)
anot <- cbind(anot, code_anot)
colnames(anot)[12] <- "code"

# Merge del archivo de anotación y de la tabla de genes

res <- merge(anot, table_genes, by.x = "code", by.y = "code")
head(res)
