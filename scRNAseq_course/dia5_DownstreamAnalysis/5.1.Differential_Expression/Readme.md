# DIA 5. Differential Gene Expression Analysis in scRNA-seq

**Servicio de Análisis Biocomputacional (SABio); CBM**\
**Edition**: May, 2026\
**Last update**: 11/05/2026

------------------------------------------------------------------------

## Introduction to differential expression analysis of scRNAseq results

Differential gene expression (DGE) analysis aims to identify genes that show significant expression differences between cell populations. In single-cell RNA sequencing (scRNA-seq), DGE analysis is commonly used to:

-   Identify marker genes for specific cell types
-   Compare transcriptional profiles between populations
-   Characterize cellular functions and biological states
-   Detect activated or disease-associated pathways

By comparing gene expression across cells or clusters, it is possible to identify genes specifically enriched in certain populations and better understand their biological roles. In this practical session, we will use **Seurat** to perform several types of differential expression analyses.

------------------------------------------------------------------------

## 0. Load libraries and data

Several R packages are required for the analysis:

-   Seurat → single-cell analysis and differential expression
-   openxlsx → export results to Excel files
-   dplyr → data manipulation and filtering
-   ggplot2 → data visualization

``` r
library(Seurat)
library(openxlsx)
library(dplyr)
library(ggplot2)

#Set random seed
set.seed(1234)

#Load annotated Seurat object
seu <- readRDS("seurat_annotated_auto.rds"

#Inspect cell type distribution
table(Idents(seu))
```

## 1. Compare one cell type versus all other cells

One common application of differential expression analysis is identifying genes specifically enriched in a given cell type. These genes are often referred to as marker genes or cell type-specific genes. As an example, we are going to calculate the marker genes of naive CD4 T cells.

In **Seurat**, the function *FindMarkers()* compares one cell population against all remaining cells. The output contains:

-   log2 fold change
-   p-values and adjusted p-values
-   percentage of expressing cells

In this example, genes with high positive log2 fold change are more highly expressed in Naive CD4 T cells compared to other populations.

``` r
#Identify marker genes for Naive CD4 T cells
deg_cd4T <- FindMarkers(
  seu,
  ident.1 = "Naive CD4 T", #Target cell population
  only.pos = TRUE #Retain only positively enriched genes
)

head(deg_cd4T)

#Export results to Excel
write.xlsx(deg_cd4T,  file = "DEGs_CD4T.xlsx", overwrite = TRUE)

#Visualize marker genes with a heatmap
DoHeatmap(
  seu,
  features = row.names(deg_cd4T)[1:30] #Plot only the first 30 markers
)
```

## 2. Compare two specific cell populations

Differential expression can also be used to directly compare two selected cell types. This analysis highlights genes that distinguish one population from another. As an example, we are going to compare CD8 T cells versus Naive CD4 T cells.

The log2 fold change for each gene is calculated relative to the reference population (Ratio between Case/Control). In this example, we use the Naive CD4 T cells population as reference, therefore:

-   log2FC \> 0 → gene overexpressed in CD8 T cells
-   log2FC \< 0 → gene underexpressed in CD8 T cells

``` r

deg_CD4vsCD8 <- FindMarkers(
  seu,
  ident.1 = "CD8 T", #Numerator group (aka "Case")
  ident.2 = "Naive CD4 T" #Denominator group (aka "control")
)

head(deg_CD4vsCD8)

#Export comparison results
write.xlsx(
  deg_CD4vsCD8,
  file = "DEGs_CD4vsCD8.xlsx",
  overwrite = TRUE
)

#Generate heatmap

DoHeatmap(
  seu,
  features = row.names(deg_CD4vsCD8)[1:30]  #Plot only the first 30 markers
)
```

## 3. Identify marker genes for all cell populations

Instead of analyzing one comparison at a time, Seurat can also identify marker genes for all cell populations simultaneously. This is useful for cluster annotation and the biological interpretation of the clustering results.

The function *FindAllMarkers()* performs differential expression analysis for every cluster or cell type. For each population the target group is compared against all remaining cells and then positively enriched genes are identified.

``` r
#Find marker genes for all cell types
degs_all <- FindAllMarkers(
  seu,
  only.pos = TRUE
)

head(degs_all)

#Export all marker genes
write.xlsx(
  degs_all,
  file = "DEGs_All_celltypes.xlsx",
  overwrite = TRUE
)

#Select top marker genes per cluster

top10 <- degs_all %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)

#Visualize the top 10 markers with a heatmap
DoHeatmap(
  seu,
  features = top10$gene
)
```
