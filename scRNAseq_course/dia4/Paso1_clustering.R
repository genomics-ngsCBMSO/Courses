#SCRIPT PARA REALIZAR EL CLUSTERING EN SEURAT

###############################
#0. CARGAR LIBRERIAS Y SEMILLA
###############################

library(Seurat)
library(scater)
library(openxlsx)
library(dplyr)

set.seed(1234)

#########################################
#1. Clustering
#########################################

seu <- readRDS("./Clustering/seurat_clustering.rds")

#K-nearest neighbour (KNN) using the same dimensions 

seu <- FindNeighbors(seu, dims = 1:13)

#Resolution: Controls the granularity of the clustering
#  - Higher values (e.g., 0.8–1.2) -> more, smaller clusters
#  - Lower values (e.g., 0.4–0.6) -> fewer, larger clusters

#By default is 0.8

seu <- FindClusters(seu, resolution=0.4)

#

jpeg("UMAP_Clusters.jpeg", units = "in", height = 10, width = 12, res = 300)
DimPlot(seu, pt.size = 1.5)
dev.off()


saveRDS(seu, "seurat_clustering.rds")

