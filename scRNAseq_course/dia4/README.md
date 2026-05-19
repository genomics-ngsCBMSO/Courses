# DIA 4. Clustering and annotation

**Servicio de Análisis Biocomputacional (SABio); CBM** \
**Edition**: May, 2026 \
**Last update**: 19/05/2026

## 🧬  Single-cell RNA-seq analysis pipeline (Day 4)

This repository contains the practical workflow to perform clustering, manual and automated annotation departing from an already dimensionality reduced scRNAseq **Seurat object**.

## 0. Prepare working session

First of all, we load libraries, set the seed (`set.seed`) to stabilize permutations and load the previously generated seurat object (`seurat_umap`)

```R
#########################################
#0. LOAD LIBRARIES AND SET SEED
#########################################

library(Seurat)
library(scater)
library(openxlsx)
library(dplyr)

set.seed(1234)

seu <- readRDS("seurat_umap.rds")
```
---

## 1. Clustering

Then, we look for KNN (K-nearest neighbour) with the same dimensions than UMAP, and, with them, look for clusters.
Clusters are visualized with DimPlot. Final Seurat object is saved.

```R
#########################################
#1. Clustering
#########################################

#K-nearest neighbour (KNN) using the same dimensions 

seu <- FindNeighbors(seu, dims = 1:13)

#Resolution: Controls the granularity of the clustering
#  - Higher values (e.g., 0.8–1.2) -> more, smaller clusters
#  - Lower values (e.g., 0.4–0.6) -> fewer, larger clusters

#By default is 0.8

seu <- FindClusters(seu, resolution=0.4)

# Save visualization of UMAP coloured by cluster

jpeg("UMAP_Clusters.jpeg", units = "in", height = 10, width = 12, res = 300)
DimPlot(seu, pt.size = 1.5)
dev.off()
```
---

## 2. Manual annotation at **cluster level** with marker genes

First, we will use `FindAllMarkers` to search positively enriched genes in each cluster (agains all other clusters).
These genes will be the cluster markers. We save the information as excel spreadsheet and visualize the top10 genes from each cluster as heatmap.

```R
#########################################
#2.1. Cluster Gene Markers
#########################################

#FindAllMarkers -> Find markers of one cluster against the rest of the cells in the dataset (only.pos = T for only positive markers)
#FindMarkers <- Find DEGs between 2 groups of cells (e.g., ident.1 = 1 / ident.2 = 2, DEGs between cluster 1 and 2)

markers <- FindAllMarkers(seu, only.pos = T)
#head(markers)

markers_list <- split(markers, markers$cluster) 
write.xlsx(markers_list, "Cluster_Markers.xlsx", overwrite = T)

#Select top10 markers in each cluster for visualization
top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)

# Draw heatmap
jpeg("Heatmap_marker.jpeg", units = "in", height = 10, width = 12, res = 300)
DoHeatmap(seu, features = top10$gene)
dev.off()
```

Manual annotation is done based on known genes and their expression patterns along clusters. 

Our practice dataset was generated from PBMCs. Already known marker genes and their corresponding cell types are:

| Marker | Celltype |
| ---   | --- |
| CD3D  | T cells |
| MS4A1 | B cells |
| LYZ   | monocytes |
| NKG7  | NK |
| PPBP  | platelets |

We plot those markers in the UMAP and peroform a comparative visualization by cluster with DotPlot, which allows us checking the percentage of cells expressing each marker gene and the average expression levels

```R
#########################################
#2.2 Visualize known gene markers
#########################################

# Coloured UMAP
jpeg("Markers.jpeg", units = "in", height = 12, width = 15, res = 300) 
FeaturePlot(
  seu,
  features = c("CD3D", "MS4A1", "LYZ", "NKG7", "PPBP"),
  cols = c("lightgrey", "blue"),
  reduction = "umap"
)
dev.off()

# Dotplot
jpeg("Dotplot_markers.jpeg", units = "in", height = 10, width = 12, res = 300)
DotPlot(
  seu,
  features = c("CD3D", "MS4A1", "LYZ", "NKG7", "PPBP")
) + RotatedAxis()
dev.off()
``` 

At this point, a biological insight is needed in order to assign identities to each cluster. 

|Cluster| Marker | Celltype |
| --- | ---   | --- |
| 0 | LYZ  | Monocytes |
| 1 | CD3D | T cells |
| 2 | CD3D | T cells |
| 3 | MS4A1 | B cells |
| 4 | LYZ  | Monocytes |
| 5 | CD3D + NKG7 | cyt T cells |
| 6 | CD3D + NKG7 | cyt T cells |
| 7 | MS4A1 | B cells |
| 8 | NKG7  | NK |
| 9 | PPBP  | Platelets |

This manual assignation is then transfered to the Seurat object. We then visualize UMAP coloured by the final manual annotation.

```R
#########################################
#2.3 Manual cluster annotation
#########################################

# Create identity vector with cell types
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

# Rename clusters and save annotation in metadata layer
seu <- RenameIdents(seu, new.cluster.ids)
seu$celltype_manual <- Idents(seu)

# Visualize annotation in UMAP
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
``` 
---

## 3. Automatic annotation at **cell level** with reference

Automatic annotation in Seurat requires a reference dataset with known cell labels. This reference can be personalized or can be downloaded from different repositories.

To load our own reference:
```R
ref <- readRDS("pbmc3k_reference.rds")
```

Seurat offers a set of references that can be consulted. 
```R
AvailableData()
```

We are interested in the PBMC reference dataset. We install and load the data:
```R
#########################################
#3.1 Reference dataset
#########################################
#Install
InstallData("pbmc3k")

#Load
reference <- LoadData("pbmc3k")
```

Both datasets have to be compatible. In order to adequately normalize the reference dataset, we use SCTransform:

```R
#########################################
#3.2 Preprocessing reference dataset
#########################################

reference <- SCTransform(reference)
```

Now we can perform automatic annotation, which is based in the search of transfer anchors. Cells from query and reference will be projected to a common PCA and KNN from both datasets will be found. Anchors link pairs of cells similar in both datasets and assign metrics to allow for identity transfer.

Then, labels are transferred from reference to query. This transfer will have a prediction score, which measures the transfer confidence based on anchors and KNNs.

```R
#########################################
#3.3 Automatic annotation
#########################################

#Find anchors: link similar cells from both datasets

anchors <- FindTransferAnchors(
  reference = reference,
  query = seu,
  normalization.method = "SCT",
  dims = 1:30
)

#Transfer reference label to each cell from query
#Pay attention to the annotation column name: seurat_annotations

predictions <- TransferData(
  anchorset = anchors,
  refdata = reference$seurat_annotations, #Nombre de la columna de la referencia donde estan los tipos celulares
  dims = 1:30
)

#Add annotation to Seurat object

seu <- AddMetaData(seu, metadata = predictions)

#Check transfer: usually "predicted.id" and "prediction.score.max"
head(seu)
```

Finally, automatic annotation can be visualized. We will colour the original UMAP according to the new annotations and to the confidence of the automatic assignation prediction.

The final object containing the cluster information, manual cluster annotation by marker genes and automatic cell annotation are saved in a new file.

```R
#########################################
#3.4 Visualization
#########################################

# Cell annotation in UMAP
jpeg("UMAP_Annotated_Automatic.jpeg", units = "in", height = 10, width = 12, res = 300)
DimPlot(
  seu,
  reduction = "umap",
  group.by = "predicted.id",
  label = TRUE,
  pt.size = 1.5
)
dev.off()

# Annotation confidence in UAMP
jpeg("Prediction_Score.jpeg", units = "in", height = 10, width = 12, res = 300)
FeaturePlot(
  seu,
  features = "prediction.score.max"
)
dev.off() 

#Save seurat object
saveRDS(seu, "seurat_annotated_auto.rds")
```







