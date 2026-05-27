# 🧬 Single-cell RNA-seq analysis (Day 3)

This repository contains the practical workflows for doublet detection, low-quality cell filtering, normalization, and dimensionality reduction.

## Doublet Detection

### Load libraries and set seed

```R
library(Seurat)
library(scDblFinder)

set.seed(1234) # ensures reproducibility of results
```

### Import CellRanger data

```R
# 1. Read the generated filtered MEX directories
data_dir <- "cellranger/PBMC/outs/filtered_feature_bc_matrix"
data <- Read10X(data.dir = data_dir)

# 2. Initialize the basic Seurat object
seu <- CreateSeuratObject(counts = data)

# Visualize the structure of the created object
head(seu)
```

### Doublet identification

```R
# 1. Convert to SingleCellExperiment format
sce <- as.SingleCellExperiment(seu)

# 2. Detect doublets
sce <- scDblFinder(sce)

# 3.  Inspect key outputs
sce$scDblFinder.class # labels each cell as singlet or doublet
sce$scDblFinder.score # probability (0–1) that a cell is a doublet

# 4. Add results to Seurat object
seu$doublet <- sce$scDblFinder.class
seu$doublet_score <- sce$scDblFinder.score

seu
head(seu)
```

### Visualization

```R
# Scatter plot
jpeg("ScatterPlot_doublets.jpeg", units = "in", height = 10, width = 10, res = 300)
FeatureScatter(seu, feature1= "nCount_RNA", feature2 = "nFeature_RNA", group.by= "doublet", cols = c("darkred", "darkblue"))
dev.off()

# Violin Plot
jpeg("ViolinPlot_doublets.jpeg", units = "in", height = 10, width = 10, res = 300)
VlnPlot(seu, features= c("nCount_RNA", "nFeature_RNA"), group.by= "doublet", cols = c("darkred", "darkblue"))
dev.off()
```

### Doublet filtering
```R
seu <- subset(seu, doublet == "singlet")
```

### Save results
```R
saveRDS(seu, "seurat_nodoub.rds")
```
