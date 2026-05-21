# DIA 5. Cell-Cell Communication Analysis

**Servicio de Análisis Biocomputacional (SABio); CBM**\
**Edition**: May, 2026\
**Last update**: 11/05/2026

------------------------------------------------------------------------

## Introduction to cell-cell communication analysis

Cells within a tissue constantly communicate with each other through signaling molecules such as cytokines, chemokines, growth factors, and membrane-bound receptors. Understanding these communication networks is essential for studying tissue organization, immune responses, development, and disease mechanisms.

Single-cell RNA sequencing (scRNA-seq) provides gene expression profiles at cellular resolution, making it possible to infer potential interactions between different cell populations based on the expression of ligands and receptors.

In this practical session, we will use **CellChat**, an R package designed to infer, analyze, and visualize intercellular communication networks from scRNA-seq data. CellChat uses curated ligand-receptor interaction databases and statistical models to: - identify potential communication events - infer signaling pathways - quantify interaction strength - visualize communication networks between cell types

## 0. Load libraries

Several R packages are required for the analysis:

-   Seurat → manipulation and analysis of single-cell datasets
-   openxlsx → export results to Excel files
-   dplyr → data manipulation and filtering
-   CellChat → inference of cell-cell communication networks
-   patchwork → combine multiple plots into complex figures

``` r
library(Seurat)
library(openxlsx)
library(dplyr)
library(CellChat)
library(patchwork)

#Set random seed 
set.seed(1234)

#Load Seurat object 
seu <- readRDS("seurat_annotated_auto.rds")
```

## 1. Create the CellChat object

Before inferring communication networks, CellChat requires as input a gene expression matrix and the metadata describing cell identities. Both can be extracted from the Seurat object. The normalized data layer is used because CellChat analyzes relative expression levels across cell populations. Cells sharing the same annotation will be grouped together during communication inference.

``` r
#Extract expression matrix 
input <- GetAssayData(seu, layer = "data")

#Create metadata table
meta <- data.frame(celltype = seu$predicted.id)
rownames(meta) <- colnames(seu)

#Inspect both tables
head(input)
head(meta)

#Create CellChat object
cellchat <- createCellChat(
  object = input, #Expression matrix
  meta = meta, #Metadata of the cells
  group.by = "celltype" #mMetadata column used to define cell populations
)

#Inspect CellChat object
cellchat
```

## 2. Load ligand-receptor interaction database

CellChat uses curated databases of known ligand-receptor interactions. These databases include: signaling pathways, ligand-receptor pairs and interaction categories. Two main databases are available:

-   CellChatDB.human → human datasets
-   CellChatDB.mouse → mouse datasets

``` r
#Select database 
CellChatDB <- CellChatDB.human 

#Explore database categories 
showDatabaseCategory(CellChatDB)

#Inspect ligand-receptor interactions
dplyr::glimpse(CellChatDB\$interaction)

#Assign database to CellChat object
cellchat@DB <- CellChatDB
```

## 3.Infer cell-cell communication networks

This section identifies biologically relevant communication events between cell populations.

### 3.1 Subset relevant genes

This step reduces the dataset to genes involved in known ligand-receptor interactions. Filtering unnecessary genes reduces computational time, improves efficiency and focuses the analysis on signaling genes.

``` r
cellchat <- subsetData(cellchat) 
```

### 3.2 Identify overexpressed genes and overexpressed ligand-receptor

CellChat identifies genes significantly enriched in each cell population. This step detects candidate signaling molecules that may participate in intercellular communication. Likewise, potential communication events are inferred by combining overexpressed ligands in one cell type and overexpressed receptors in another cell type. Only biologically plausible ligand-receptor pairs are retained.

``` r
cellchat <- identifyOverExpressedGenes(cellchat)

cellchat <- identifyOverExpressedInteractions(cellchat)
```

### 3.3 Compute communication probabilities and infer signaling pathways

This step estimates the probability and strength of communication between cell populations. The probability is calculated using the expression levels, the ligand-receptor abundance and the interaction models.

Interactions are then aggregated into signaling pathways (Examples: MHC-II, TGF-beta, IFN signaling). This allows pathway-level interpretation of communication networks. Finally CellChat summarizes communication events across all cell populations, with the resulting final network containing:

-   The number of interactions

-   The strength of the identified interactions

-   Pathway-level communication patterns

``` r
#Calculate probabilities
cellchat <- computeCommunProb(cellchat)

#Optional filtering
#cellchat <- filterCommunication(cellchat, min.cells = 10) 
#This optional step removes interactions involving very small cell populations. Filtering low-cell groups may reduce noise and spurious interactions.

#Infer signaling pathways 
cellchat <- computeCommunProbPathway(cellchat) 

#Aggregate communication network 
cellchat <- aggregateNet(cellchat)
```

## 4.Visualize communication networks

### 4.1 Plot number of interactions

With the function *netVisual_circle* we generate a circular plot that shows which cell types comunnicate and how many interactions occur between them. Thicker edges represent more communication events.

``` r
netVisual_circle(cellchat@net$count, 
                 weight.scale = TRUE, 
                 label.edge = FALSE, 
                 title.name = "Number of interactions"
```

### 4.2 Plot interaction strength

In addition, the same circular plot can be made to represent the overall intensity of communication between cell populations.

``` r
netVisual_circle(cellchat@net$weight, 
                  weight.scale = TRUE,
                  label.edge = FALSE,
                  title.name = "Interaction strength")
```

## 5.Visualize specific signaling pathways

Specific signaling pathways can be explored individually. As an example, we visualize MHC-II signaling network, which is involved in antigen presentation and immune activation.

``` r
#Inspect detected pathways
#This displays the signaling pathways identified during the analysis.

head(cellchat@netP$pathways) 

#Visualize MHC-II signaling network 
#Circle plot 
netVisual_aggregate( cellchat, 
                     signaling = "MHC-II", 
                     layout = "circle")
                     
#Chord diagram 
netVisual_aggregate( cellchat,
                     signaling = "MHC-II", 
                     layout = "chord")

#Heatmap visualization 
netVisual_heatmap(cellchat, 
                  signaling = "MHC-II", 
                  color.heatmap = "Reds")
```

## 6. Save results Saving the object preserves the inferred interactions,

pathways and networks, allowing the analysis to be reloaded without recomputing all steps.

``` r
saveRDS(cellchat,
    "cellchat_object.rds") 
```
