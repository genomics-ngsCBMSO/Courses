# DIA 5. Trajectory Analysis in scRNA-seq

**Servicio de Análisis Biocomputacional (SABio); CBM**\
**Edition**: April, 2026\
**Last update**: 11/05/2026

------------------------------------------------------------------------

## Introduction to trajectory analysis

In many biological systems, cells do not exist as static populations. Instead, they transition through continuous biological processes such as cell differentiation, immune activation or disease progression.

Trajectory analysis aims to reconstruct these dynamic transitions using single-cell RNA sequencing (scRNA-seq) data. Rather than grouping cells into isolated clusters, trajectory methods attempt to organize cells along continuous developmental paths and infer transitional cellular states, in order to study the progression through a biological process.

One of the most widely used approaches is **pseudotime analysis**. Pseudotime is a computational measure that estimates the relative progression of cells through a biological trajectory.

Important!:

-   pseudotime is **not real chronological time**

-   it represents the inferred order of cellular transitions based on gene expression similarity

Cells with similar transcriptional profiles are placed close together, while progressively changing cells are arranged along a trajectory.

## 0. Load libraries and data

Several R packages are required for the analysis:

-   openxlsx → export results to Excel files
-   dplyr → data manipulation
-   patchwork → combine visualizations
-   SingleCellExperiment → standardized single-cell data structure
-   monocle3 → trajectory and pseudotime inference
-   ggplot2 → plotting and visualization

``` r
library(openxlsx)
library(dplyr)
library(patchwork)
library(SingleCellExperiment)
library(monocle3)
library(ggplot2)

#Set random seed
set.seed(1234)

#Load trajectory dataset
cds <- readRDS("trayectory_dataset_CDS.rds")

#Inspect dataset structure
cds
rowData(cds) #metadata associated with genes
colData(cds) #metadata associated with cells
```

## 0PTIONAL: EXAMPLE SEURAT TO MONOCLE3

If we want to work with an already analyzed Seurat dataset, we can import the clusters and reductions to the monocle object

```r
# 1. Get data from Seurat & create el CDS object
expression_matrix <- GetAssayData(seu, assay = "RNA", slot = "counts")
cell_metadata <- seu@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix), 
                              row.names = rownames(expression_matrix))

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

# 2. Transfer UMAP to Monocle 3
reducedDims(cds)[["UMAP"]] <- Embeddings(seu, reduction = "umap")

# 3. Transfer Clusters to Monocle3
cds@clusters[["UMAP"]][["clusters"]] <- Idents(seu)

# 4. Create Partitions
names_cells <- rownames(seu@meta.data)
partitions <- factor(rep(1, length(names_cells)), levels = "1")
names(partitions) <- names_cells
cds@clusters[["UMAP"]][["partitions"]] <- partitions

#We can continue here straight to the Section 2. Trajectory inference
```

## 1. Dimensionality reduction and clustering

Trajectory analysis begins by organizing cells according to transcriptional similarity.

``` r
#Reduce dimensionality with UMAP and visualize
cds <- reduce_dimension(cds)

plot_cells(
  cds,
  cell_size = 0.5
) +
theme(
  axis.text = element_text(size = 15),
  axis.title = element_text(size = 15)
)

#Visualize cells by cell type
plot_cells(
  cds,
  cell_size = 0.5,
  color_cells_by = "cell.type",
  group_label_size = 5
) +
theme(
  axis.text = element_text(size = 15),
  axis.title = element_text(size = 15)
)

#Cluster cells
cds <- cluster_cells(cds)

plot_cells(
  cds,
  cell_size = 0.5,
  color_cells_by = "cluster",
  group_label_size = 7
) +
theme(
  axis.text = element_text(size = 15),
  axis.title = element_text(size = 15)
)
```

## 2. Trajectory inference

### 2.1 Construct trayectory graph

The next step is reconstructing the biological trajectories connecting cells.

Monocle3 constructs a graph representing potential cellular transitions, which forms the basis for pseudotime analysis. The graph connects cells according to transcriptional similarity and attempts to infer developmental paths, branching events and lineage relationships among the cells.

``` r
#Learn trajectory graph
cds <- learn_graph(cds)

#Visualize trajectory structure
plot_cells(
  cds,
  color_cells_by = "cell.type",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  group_label_size = 3
)
```

### 2.2. Define the trajectory root

Pseudotime always requires defining a biological starting point. This starting point represents the earliest or most initial cellular state in the process being studied. Visualization of the trajectory structure helps identify the most appropriate starting population.

One option is selecting a known biological population as the root, for example progenitor populations or a population of naive immune cells. Pseudotime is calculated relative to those cells

``` r
root_cells <- colnames(cds)[colData(cds)$celltype == "Monocytes"]
cds <- order_cells(cds, root_cells = root_cells)
```

However, Monocle3 also allows to select the root interactively, which allows to select manually one specific node of the trajectory graph as the root node

``` r
cds <- order_cells(cds)

#Visualize pseudotime
plot_cells(
  cds,
  cell_size = 0.5,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)
```

## 3. Identify genes changing along pseudotime

One major goal of trajectory analysis is identifying genes whose expression changes dynamically during the biological process. Dynamic genes help interpret the molecular programs underlying cellular transitions.

In Monocle3, the function *graph_test()* identifies genes whose expression is associated with the trajectory structure. The method evaluates whether gene expression varies significantly across neighboring cells along the graph. The output contains: statistical significance values, q-values and trajectory-associated expression metrics.

Genes with low q-values are strongly associated with pseudotime progression.

``` r
#Detect dynamic genes (This step may take some time depending on dataset size)
deg_pseudotime <- graph_test(
  cds,
  neighbor_graph = "principal_graph"
)

head(deg_pseudotime)

#Select top dynamic genes
top_genes <- rownames(
  subset(deg_pseudotime, q_value < 0.05)
)[1:6]

#Visualize dynamic genes
plot_cells(
  cds,
  genes = top_genes,
  show_trajectory_graph = TRUE
)
```

## 6. Save results

Save Monocle3 object

``` r
saveRDS(cds, "monocle3_trajectory.rds")
```
