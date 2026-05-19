#SCRIPT PARA REALIZAR LA ANOTACION MANUAL POR MARCADORES EN SEURAT

###############################
#0. CARGAR LIBRERIAS Y SEMILLA
###############################

library(Seurat)
library(scater)
library(openxlsx)
library(dplyr)

set.seed(1234)

seu <- readRDS("seurat_clustering.rds")

#########################################
#1. Cluster Gene Markers
#########################################

#FindAllMarkers -> Find markers of one cluster against the rest of the cells in the dataset (only.pos = T for only positive markers)
#FindMarkers <- Find DEGs between 2 groups of cells (e.g., ident.1 = 1 / ident.2 = 2, DEGs between cluster 1 and 2)

markers <- FindAllMarkers(seu, only.pos = T)
#head(markers)

markers_list <- split(markers, markers$cluster) 
write.xlsx(markers_list, "Cluster_Markers.xlsx", overwrite = T)

#Seleccionamos los top10 marcadores de cada cluster para plotearlos
top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)

jpeg("Heatmap_marker.jpeg", units = "in", height = 10, width = 12, res = 300)
DoHeatmap(seu, features = top10$gene)
dev.off()


#########################################
#2. Anotación manual por genes marcadores
#########################################

#2.1 Visualización de genes marcadores conocidos

#En scRNA-seq, la anotacion manual se basa en:
# - genes conocidos (marcadores)
# - su patron de expresion en los clusters


#Algunos marcadores clasicos en PBMC:
#CD3D  -> celulas T
#MS4A1 -> celulas B
#LYZ   -> monocitos
#NKG7  -> NK
#PPBP  -> plaquetas

jpeg("Markers.jpeg", units = "in", height = 12, width = 15, res = 300) 
FeaturePlot(
  seu,
  features = c("CD3D", "MS4A1", "LYZ", "NKG7", "PPBP"),
  cols = c("lightgrey", "blue"),
  reduction = "umap"
)
dev.off()

#2.2 Visualizacion comparativa por cluster

#DotPlot permite ver:
# - porcentaje de celulas que expresan el gen
# - nivel medio de expresion

jpeg("Dotplot_markers.jpeg", units = "in", height = 10, width = 12, res = 300)
DotPlot(
  seu,
  features = c("CD3D", "MS4A1", "LYZ", "NKG7", "PPBP")
) + RotatedAxis()
dev.off()


#2.3 Interpretacion biologica


#En este punto del analisis, interpretamos manualmente:

#Cluster	Marker	        Anotación Manual
#0      	LYZ   	          Monocytes
#1      	CD3D   	          T cells
#2      	CD3D   	          T cells		
#3		    MS4A1  	          B cells
#4      	LYZ    	          Monocytes
#5		    CD3D + NKG7	      cyt T cells
#6		    CD3D + NKG7	      cyt T	cells
#7		    MS4A1  	          B cells		
#8		    NKG7    	        NK
#9        PPBP   	          Platelets

#2.4 Asignacion manual de identidades


#Creamos un vector con los nombres de los tipos celulares
#El orden DEBE coincidir con levels(seu)

levels(seu)

new.cluster.ids <- c(
  "Monocytes",
  "T cells",
  "T cells",
  "B cells",
  "Monocytes",
  "cyt T cells",
  "cyt T cells",
  "B cells",
  "NK cells",
  "Platelets"
)

names(new.cluster.ids) <- levels(seu)

#Renombramos los clusters
seu <- RenameIdents(seu, new.cluster.ids)

#Guardamos la anotacion en metadata
seu$celltype_manual <- Idents(seu)


#1.5 Visualizacion final anotada

jpeg("UMAP_Annotated_Manual.jpeg", units = "in", height = 10, width = 12, res = 300)

DimPlot(
  seu,
  reduction = "umap",
  group.by = "celltype_manual",
  label = TRUE,
  label.box = TRUE,
  pt.size = 1.5
)

dev.off()

#########################################
# Guardar objeto anotado
#########################################

saveRDS(seu, "seurat_markeranot.rds")





