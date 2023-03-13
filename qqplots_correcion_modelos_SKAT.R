rm(list=ls())
setwd("C:/Users/USUARIO/Desktop/FPIES")
library(Matrix)
library(SPAtest)
library(RSpectra)
library(SKAT)
library(qqman)

# Lectura de archivos de la cohorte FPIES

skat <- read.table("file_model_skat.txt", header = TRUE)
burden <- read.table("file_model_burden.txt", header = TRUE)
skato <- read.table("file_model_skato.txt", header = TRUE)

# Lectura de archivos de la cohorte FPIES-A

skat_a <- read.table("file_model_skat.txt", header = TRUE)
burden_a <- read.table("file_model_burden.txt", header = TRUE)
skato_a <- read.table("file_model_skato.txt", header = TRUE)

# Lectura de archivos de la cohorte FPIES-B

skat_b <- read.table("file_model_skat.txt", header = TRUE)
burden_b <- read.table("file_model_burden.txt", header = TRUE)
skato_b <- read.table("file_model_skato.txt", header = TRUE)

# QQ-plots 

qq(skat$P.value, main = "QQ PLOT SKAT - PONDERADO FPIES")
qq(skato$P.value, main = "QQ PLOT SKATO - PONDERADO FPIES ")
qq(burden$P.value, main = "QQ PLOT BURDEN - PONDERADO FPIES")

qq(skat_a$P.value, main = "QQ PLOT SKAT - PONDERADO FPIES-A")
qq(skato_a$P.value, main = "QQ PLOT SKATO - PONDERADO FPIES-A")
qq(burden_a$P.value, main = "QQ PLOT BURDEN - PONDERADO FPIES-A")

qq(skat_b$P.value, main = "QQ PLOT SKAT - PONDERADO FPIES-B")
qq(skato_b$P.value, main = "QQ PLOT SKATO - PONDERADO FPIES-B")
qq(burden_b$P.value, main = "QQ PLOT BURDEN - PONDERADO FPIES-B")

# La corrección de Bonferroni se hace siempre entre el número de comparaciones que
# tengas, es decir, en estos modelos sería 0.05/num genes. El número de genes es el
# número de filas del archivo de resultados

# Corrección de Bonferroni FPIES vs.IBS

alpha_adj <- 0.05/4526

skat_adj <- skat[which(skat$P.value<alpha_adj),]
burden_adj <- burden[which(burden$P.value<alpha_adj),]
skato_adj <- skato[which(skato$P.value<alpha_adj),]

# Corrección de Bonferroni FPIES-A vs.IBS

alpha_adj <- 0.05/3387

skat_hc_adj <- skat_hc[which(skat_hc$P.value<alpha_adj),]
burden_hc_adj <- burden_hc[which(burden_hc$P.value<alpha_adj),]
skato_hc_adj <- skato_hc[which(skato_hc$P.value<alpha_adj),]

# Corrección de Bonferroni FPIES-B vs.IBS

alpha_adj <- 0.05/3530

skat_lc_adj <- skat_lc[which(skat_lc$P.value<alpha_adj),]
burden_lc_adj <- burden_lc[which(burden_lc$P.value<alpha_adj),]
skato_lc_adj <- skato_lc[which(skato_lc$P.value<alpha_adj),]

