# 🧬 Single-cell RNA-seq analysis (Day 3)

This repository contains the practical workflows for doublet detection, low-quality cell filtering, normalization, and dimensionality reduction.

## 🔍  Doublet Detection

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

## 🧹 Low-quality cell filtering

### Load libraries and set seed

```R
library(Seurat)
library(scater)

set.seed(1234) # ensures reproducibility of results
```

### Load Seurat object (after doublet removal)

```R
seu <- readRDS("seurat_nodoub.rds")
```

### Quality control metrics

```R
# Inspect metadata
head(seu)

# Calculate mitochondrial gene percentage
sum(grepl("^MT-", rownames(seu)))  # how many genes match the MT- pattern
seu$percent.mt <- PercentageFeatureSet(seu, pattern = "^MT-")
```

### Visualize QC metrics

```R
jpeg("ViolinPlot_quality.jpeg", units = "in", height = 10, width = 10, res = 300)
VlnPlot(seu, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
dev.off()
```

### Log-transformed metrics

```R
seu$log10_nCount <- log10(seu$nCount_RNA)
seu$log10_nFeature <- log10(seu$nFeature_RNA)

jpeg("ViolinPlot_quality_log10.jpeg", units = "in", height = 10, width = 10, res = 300)
VlnPlot(seu, features = c("log10_nCount", "log10_nFeature", "percent.mt"))
dev.off()
```

### Detect outliers

```R
seu$nCount_outliers <- isOutlier(seu$nCount_RNA, nmads = 5)
seu$nFeature_outliers <- isOutlier(seu$nFeature_RNA, nmads = 5)
seu$mt_outliers <- isOutlier(seu$percent.mt, nmads = 5, type = "higher")

seu$outliers <- "KEEP"
seu$outliers[seu$nCount_outliers == TRUE |
             seu$nFeature_outliers == TRUE |
             seu$mt_outliers == TRUE] <- "FILTER" # Cells outside normal ranges are marked as FILTER
```


### Visualize outliers

```R
p1 <- FeatureScatter(seu, group.by = "outliers",
                     feature1 = "nCount_RNA", feature2 = "percent.mt",
                     cols = c("darkred", "darkblue"))

p2 <- FeatureScatter(seu, group.by = "outliers",
                     feature1 = "nFeature_RNA", feature2 = "percent.mt",
                     cols = c("darkred", "darkblue"))

p3 <- FeatureScatter(seu, group.by = "outliers",
                     feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                     cols = c("darkred", "darkblue"))

jpeg("Scatter_outliers.jpeg", units = "in", height = 12, width = 15, res = 300)
p1 + p2 + p3
dev.off()
```


### Alternative visualization (ggplot2)

```R
library(ggplot2)
library(tidyr)
library(dplyr)

# Extract data
df <- FetchData(seu, vars = c("nCount_RNA", "nFeature_RNA", "percent.mt", "outliers"))

# Convert to long format
df_long <- df %>%
  pivot_longer(cols = c(nCount_RNA, nFeature_RNA, percent.mt),
               names_to = "feature",
               values_to = "value")

# Plot
ggplot(df_long, aes(x = feature, y = value)) +
  geom_violin(fill = "lightgrey") +
  geom_jitter(aes(color = outliers), width = 0.2, size = 1) +
  facet_wrap(~feature, scales = "free") +
  theme_classic() +
  scale_color_manual(values = c("KEEP" = "darkblue", "FILTER" = "darkred"))
```

### Filter low-quality cells
```R
seu <- subset(seu, outliers == "KEEP")
```

### Save results
```R
saveRDS(seu, "seurat_filtered.rds")
```
