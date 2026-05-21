# DIA 5. Cell Type Proportion Analysis in scRNA-seq

**Servicio de Análisis Biocomputacional (SABio); CBM**\
**Edition**: May, 2026\
**Last update**: 11/05/2026

------------------------------------------------------------------------

## Introduction to cell type proportion analysis

One important application of single-cell RNA sequencing is the comparison of cellular composition between biological conditions or samples. Cell type proportion analysis evaluates whether the relative abundance of cell populations differs significantly between samples or experimental groups. Changes in the abundance of specific cell populations may reflect immune activation, treatment response or developmental changes for example.

Unlike previous analyses performed on one sample, proportion analysis requires comparison across samples or conditions. For this reason, on this practical session, we need to use a different Seurat object that contains cells from multiple donors combined into a single dataset.

In this example:

-   Donor1 → first biological sample
-   Donor2 → second biological sample

The object already contains the normalized expression data of both donors combinded, cell annotations and the donor/sample metadata.

## 0. Load libraries and data

Several R packages are required for the analysis:

-   Seurat → handling scRNA-seq datasets
-   openxlsx → export results to Excel files
-   dplyr → data manipulation
-   ggplot2 → data visualization

``` r
library(Seurat)
library(openxlsx)
library(dplyr)
library(ggplot2)

#Set random seed
set.seed(1234)

#Load combined Seurat object
seu <- readRDS("combined_donors.rds")

head(seu)

#Count cells per donor
table(seu$Sample)
```

## 1. Extract metadata for proportion analysis

To compare cellular composition between donors, we need the cell type annotations (predicted cell identity) and the sample/donor origin information These variables are stored in the Seurat metadata.

``` r
#Extract relevant metadata columns
md <- seu@meta.data %>%
      select(CellType = predicted.id, Sample = Sample)


CellType → predicted cell identity
Sample → donor/sample origin

#Convert variables to factors
md$CellType <- factor(md$CellType)
md$Sample <- factor(md$Sample)

head(md)
```

## 2. Global comparison of cell type frequencies

The first analysis evaluates whether overall cellular composition differs between donors.

For this, we create a contingency table that summarizes the abundance of each cell type in each donor. Then we can calculate whether the distribution of cell types differs significantly between samples with a chi-square test:

-   Null hypothesis: cell type proportions are similar across donors
-   Alternative hypothesis: at least one cell population differs in abundance

A p-value belowe 0.05 indicates statistically significant differences in cell type proportions between the two donors.

An important aspect to take into account is that in scRNA-seq data, some cell populations may contain very few cells. Low expected frequencies can violate the assumptions of the classical chi-square test.

Within the function *chisq.test* in base R we use the parameter *simulate.p.value = TRUE*, which computes the p-value using Monte Carlo simulations instead of relying on theoretical approximations. This produces more robust p-values and a better performance for sparse contingency tables.

``` r
#Contingency table
freqs <- table(md$Sample, md$CellType)
freqs

#Perform chi-square test
chisq_res <- chisq.test(freqs, simulate.p.value=TRUE)
chisq_res
```

## 3. Cell type-specific proportion tests

The global chi-square test indicates whether differences exist overall, but it does not specify which cell types are significantly altered. To identify specific populations with differential abundance, individual statistical tests are performed for each cell type.

For this step we load a custom R function that performs cell type-specific proportion testing. The function compares the abundance of each cell population between donors independently.

``` r
#Load custom testing function
source("ctprop_test.R")

#Convert frequency table into data frame
freqs_df <- as.data.frame(freqs)
colnames(freqs_df) <- c("Sample", "CellType", "Frequency")
Explanation

#Calculate proportions
#For each sample: cell type frequency / total number of cells

freqs_df$Proportion <- freqs_df$Frequency / ave(freqs_df$Frequency,
  freqs_df$Sample,
  FUN = sum
)

#Run cell type-specific tests
results <- ctprop_test(freqs_df, freqs)

#Export statistical results
write.xlsx(results, "DiffProp_results.xlsx", row.names = FALSE)
```

## 4. Visualization of cell type proportions

``` r
#Barplot visualization
ggplot(freqs_df, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(label = scales::percent(Proportion, accuracy = 0.1)),
    position = position_stack(vjust = 0.5),
    size = 3,
    fontface = "bold"
  ) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Proporción de tipos celulares por muestra",
    y = "Proporción",
    x = "Sample"
  ) +
  theme_minimal(base_size = 15)

#Pie chart visualization
ggplot(freqs_df, aes(x = "", y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", width = 1) +
  geom_text(
    aes(label = scales::percent(Proportion, accuracy = 0.1)),
    position = position_stack(vjust = 0.5),
    size = 3,
    fontface = "bold"
  ) +
  coord_polar(theta = "y") +
  facet_wrap(~ Sample) +
  theme_void(base_size = 15) +
  labs(title = "Composición celular por muestra")
```
