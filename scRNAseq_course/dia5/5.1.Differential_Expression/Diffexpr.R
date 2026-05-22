#SCRIPT PARA REALIZAR ANALISIS DE EXPRESIÓN DIFERENCIAL ENTRE TIPOS CELULARES

###############################
#0. CARGAR LIBRERIAS Y SEMILLA
###############################

library(Seurat)
library(openxlsx)
library(dplyr)
library(ggplot2)

set.seed(1234)

seu <- readRDS("seurat_annotated_auto.rds")

#Usamos la anotacion automatica como identidades celulares
Idents(seu) <- "predicted.id"

table(Idents(seu))

#########################################
#1. Comparar un tipo celular vs el resto
#########################################

#Ejemplo: genes especificos de monocitos

deg_cd4T <- FindMarkers(
  seu,
  ident.1 = "Naive CD4 T",
  only.pos = TRUE
)

head(deg_cd4T)

write.xlsx(deg_cd4T,  file = "DEGs_CD4T.xlsx", overwrite = T)

jpeg("Heatmap_DEGs_CD4T.jpeg", units = "in", height = 12, width = 15, res = 300)
DoHeatmap(seu, features = row.names(deg_cd4T)[1:30])
dev.off()

#########################################
#2 Comparar dos tipos celulares concretos
#########################################

#Ejemplo: T cells vs B cells

deg_CD4vsCD8<- FindMarkers(
  seu,
  ident.1 = "CD8 T", #Identidad del caso (numerador)
  ident.2 = "Naive CD4 T" #Identidad de la referencia (denominador)
)

#Si un gen tiene log2FC < 0 , está infra expresado en CD8 T
#Si un gen tiene log2FC > 0, está sobre expresado en CD8 T

head(deg_CD4vsCD8)

write.xlsx(deg_cd4T,  file = "DEGs_CD4vsCD8.xlsx", overwrite = T)

jpeg("Heatmap_DEGs_CD4vsCD8.jpeg", units = "in", height = 12, width = 15, res = 300)
DoHeatmap(seu, features = row.names(deg_CD4vsCD8)[1:30])
dev.off()

###################################################
#3. Genes marcadores de TODOS los tipos celulares
###################################################

degs_all <- FindAllMarkers(
  seu,
  only.pos = TRUE
)

head(degs_all)

write.xlsx(deg_cd4T,  file = "DEGs_All_celltypes.xlsx", overwrite = T)

top10 <- degs_all %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)

jpeg("Heatmap_DEGs_All_celltypes.jpeg", units = "in", height = 12, width = 15, res = 300)
DoHeatmap(seu, features = top10$gene)
dev.off()



