rm(list=ls())
setwd("C:/Users/USUARIO/Desktop/FPIES")

# Lectura de archivos
fam <- read.table("Total_FPIES/total_fpies_assoc/total_fpies_vs_khor.fam",header = FALSE, 
                  col.names = c("FID","IID", "V3", "V4", "V5", "V6"))
samples_khor <-read.csv("KFE/khor_sample_sheet_2.csv",sep = ";",header=TRUE)

# Diferenciar casos de controles
fpies <- fam[which(startsWith(fam[,1], "18") | startsWith(fam[,1], "19")| startsWith(fam[,1], "A") | 
                     startsWith(fam[,1], "B") ),]

fam[,6] <- ifelse(fam[,1]%in%fpies[,1], 2, 1)

# Escritura archivos
write.table(fam, "Total_FPIES/total_fpies_assoc/total_fpies_vs_khor.fam", row.names = FALSE,
            col.names = FALSE, quote = FALSE)
