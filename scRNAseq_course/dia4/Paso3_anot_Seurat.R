#SCRIPT PARA REALIZAR LA ANOTACION CELULAR AUTOMATICA CON SEURAT

###############################
#0. CARGAR LIBRERIAS Y SEMILLA
###############################

library(Seurat)
library(SeuratData)
library(scater)
library(openxlsx)
library(dplyr)

set.seed(1234)

seu <- readRDS("seurat_markeranot.rds")

#########################################
#1. DATASET DE REFERENCIA
#########################################

##QUE NECESITAMOS?

#La anotacion automatica en Seurat requiere:
# - un dataset de referencia
# - con etiquetas celulares conocidas

#Si tuvieramos un dataset descargado lo cargariamos como un dataset normal:
#ref <- readRDS("pbmc3k_reference.rds")

#Seurat ofrece una serie de referencias 

AvailableData()

#En nuestro caso nos interesa instalar la referencia para el dataset pbmc
InstallData("pbmc3k")

#Lo cargamos
reference <- LoadData("pbmc3k")

#########################################
#2. PREPROCESAMIENTO
#########################################

#Ambos datasets deben estar normalizados de forma compatible
#Lo mas recomendable: usar SCTransform en ambos

reference <- SCTransform(reference)

#########################################
#3. ANOTACION AUTOMATICA
#########################################

#Los anchors conectan celulas similares entre datasets

anchors <- FindTransferAnchors(
  reference = reference,
  query = seu,
  normalization.method = "SCT",
  dims = 1:30
)

#Se transfiere la etiqueta del dataset de referencia a cada celula del dataset query
#Es importante saber el nombre de la columna con la anotación, en este caso: seurat_annotations

predictions <- TransferData(
  anchorset = anchors,
  refdata = reference$seurat_annotations, #Nombre de la columna de la referencia donde estan los tipos celulares
  dims = 1:30
)

#Añadir anotacion al objeto Seurat

seu <- AddMetaData(seu, metadata = predictions)

#La columna por defecto suele ser "predicted.id" y "prediction.score.max" (confianza)
head(seu)

#4.6 Visualizacion

jpeg("UMAP_Annotated_Automatic.jpeg", units = "in", height = 10, width = 12, res = 300)

DimPlot(
  seu,
  reduction = "umap",
  group.by = "predicted.id",
  label = TRUE,
  pt.size = 1.5
)

dev.off()

#Evaluar confianza de la anotacion
#no todas las predicciones son igual de fiables

jpeg("Prediction_Score.jpeg", units = "in", height = 10, width = 12, res = 300)
FeaturePlot(
  seu,
  features = "prediction.score.max"
)
dev.off() 

#########################################
#Guardar objeto
#########################################

saveRDS(seu, "seurat_annotated_auto.rds")


