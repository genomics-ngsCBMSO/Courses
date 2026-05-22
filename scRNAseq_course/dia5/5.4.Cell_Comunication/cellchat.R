#SCRIPT PARA REALIZAR ANALISIS DE COMUNICACION CELULAR CON CELLCHAT

###############################
#0. CARGAR LIBRERIAS Y SEMILLA
###############################

library(Seurat)
library(openxlsx)
library(dplyr)
library(CellChat)
library(patchwork)

set.seed(1234)

seu <- readRDS("seurat_annotated_auto.rds")

#########################################
#1. CREAR EL OBJETO CELLCHAT
#########################################

#Extraer matriz de expresión
input <- GetAssayData(seu, layer = "data")

#Metadata con tipos celulares
meta <- data.frame(celltype = seu$predicted.id)
rownames(meta) <- colnames(seu)

head(input)
head(meta)

cellchat <- createCellChat(
  object = input,
  meta = meta,
  group.by = "celltype"
)

cellchat

str(cellchat)

#########################################
#2. BASE DE DATOS DE INTERACCIONES
#########################################

#CellChat usa bases de datos de ligando-receptor conocidas.
#Para humano -> CellChatDB.human
#Para ratón -> CellChatDB.mouse

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

#Asignamos la base de datos al objeto
cellchat@DB <- CellChatDB

#########################################
#3. IDENTIFICAR INTERACCIONES
#########################################

#3.1: Filtrar genes relevantes
cellchat <- subsetData(cellchat)

#3.2: Detectar genes sobreexpresados
cellchat <- identifyOverExpressedGenes(cellchat)

#3.3: Detectar interacciones ligando-receptor relevantes
cellchat <- identifyOverExpressedInteractions(cellchat)

#3.4: Calcular probabilidad de comunicación entre tipos celulares
cellchat <- computeCommunProb(cellchat)
# OPCIONAL: Filtrar los eventos de comunicación si hay un número de células bajo en alguno de los grupos
#cellchat <- filterCommunication(cellchat, min.cells = 10)

#3.5: Inferir comunicación a nivel de vías (pathways)
cellchat <- computeCommunProbPathway(cellchat)

#3.6: Agregar red de comunicación
cellchat <- aggregateNet(cellchat)

#3.7: Plotear los resultados de interacción

jpeg("Cellchat_numberinteractions.jpeg", units = "in", height = 10, width = 6, res = 300)
netVisual_circle(cellchat@net$count, weight.scale = TRUE, label.edge = FALSE,
                 title.name = "Number of interactions")
dev.off()

jpeg("Cellchat_strenghtinteractions.jpeg", units = "in", height = 10, width = 6, res = 300)
netVisual_circle(cellchat@net$weight, weight.scale = TRUE, label.edge = FALSE,
                 title.name = "Interaction strength")
dev.off()

#########################################
#4. VISUALIZAR PATHWAYS 
#########################################

#Podemos ver las rutas que se han obtenido en el analisis 
head(cellchat@netP$pathways)

#Visualizamos MHC-II (Complejo Mayor de Histocompatibilidad clase II)

jpeg("Cellchat_MHCII_circle.jpeg", units = "in", height = 10, width = 6, res = 300)
netVisual_aggregate(cellchat, signaling = "MHC-II", layout = "circle")
dev.off()

jpeg("Cellchat_MHCII_chrod.jpeg", units = "in", height = 10, width = 6, res = 300)
netVisual_aggregate(cellchat, signaling = "MHC-II", layout = "chord")
dev.off()

jpeg("Cellchat_MHCII_heatmap.jpeg", units = "in", height = 10, width = 6, res = 300)
netVisual_heatmap(cellchat, signaling = "MHC-II", color.heatmap = "Reds")
dev.off()

#########################################
#5. GUARDAR OBJETO Y SALIR  
#########################################

saveRDS(cellchat, "cellchat_object.rds")




