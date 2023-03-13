rm(list=ls())
setwd("/path/to/directory")
library(Matrix)
library(SPAtest)
library(RSpectra)
library(SKAT)

file_bed <- "file.bed"
file_bim <- "file.bim"
file_fam <- "file.fam"
file_setid <- "file.SetID"

# Primero, se crean los ficheros SSD e Info

Generate_SSD_SetID(file_bed, file_bim, file_fam, file_setid, "file.SSD",
                   "file.info", Is.FlipGenotype = TRUE)

# Después, se abre el fichero SSD, el fichero Info y el fichero de los pesos

fam <- Read_Plink_FAM(Filename = file_fam, Is.binary = FALSE)
y <- fam$Phenotype
y <- (y-1)

pesos<-Read_SNP_WeightFile("file.txt")
ssd_info <- Open_SSD(File.SSD = "file.SSD", File.Info = "file.info" )

# Número de muestras

ssd_info$nSample

# Número de sets

ssd_info$nSets

# Generar modelo y residuos para SKAT. Se pueden generar varios tipos de modelos jugando
# con los parámetros y ver cual de todos ellos se adapta mejor a los datos

obj <- SKAT_Null_Model(y ~ 1, out_type = "D")
out <- SKAT.SSD.All(SSD.INFO = ssd_info, obj) 
out
out_maf <- SKAT.SSD.All(SSD.INFO = ssd_info, obj, max_maf = 0.05)
out_pesos <- SKAT.SSD.All(SSD.INFO = ssd_info, obj, obj.SNPWeight = pesos)
out_maf_pesos <- SKAT.SSD.All(SSD.INFO = ssd_info, obj, obj.SNPWeight = pesos, max_maf = 0.05)

# El output del objeto SKAT.SSD.All tiene un output dataframe llamado "results"

output_df <- out$results

# Eliminar  todas aquellas muestras en donde la columna N.Marker.Test (numero de
# marcadores a testear para la asociación después de excluir los no polimorficos 
# o los que presentan alta tasa de missing-rate) sea < 7.

out_test <- output_df[which(output_df$N.Marker.Test>7),]

# Los resultados están ordenados por el SetID y no por el p-valor; por lo que para
# encontrar los genes que están más asociados con el fenotipo, se deben ordenar por 
# p-valor

head(out_test[order(out_test$P.value),])
out_adj <- (out_test[order(out_test$P.value),])

# se puede guardar en un fichero con la función write.table
write.table(out_adj, file = "file_save.txt", col.names = TRUE, row.names = FALSE, 
            quote = FALSE)

# Se prueban otros tipos de test: Burden y SKATO (se pueden hacer las mismas pruebas
# que para SKAT)

out_burden <- SKAT.SSD.All(SSD.INFO = ssd_info, obj, r.corr = 1 )
out_skato <- SKAT.SSD.All(SSD.INFO = ssd_info, obj, method = "SKATO")

burden_df <- out_burden$results
skato_df <- out_skato$results

burden_test <- burden_df[which(burden_df$N.Marker.Test>7),]
skato_test <- skato_df[which(skato_df$N.Marker.Test>7),]

burden_adj <- (burden_test[order(burden_test$P.value),])
head(burden_adj) 

skato_adj <- (skato_test[order(skato_test$P.value),])
head(skato_adj)

# Se guardan los resultados

write.table(burden_adj, file = "file_save.txt", col.names = TRUE, row.names = FALSE, 
            quote = FALSE)
write.table(skato_adj, file = "file_save.txt", col.names = TRUE, row.names = FALSE,
            quote = FALSE)

# Después de usar los ficheros SSD, se deben cerrar

Close_SSD()
