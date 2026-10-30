#!/usr/bin/env Rscript
## R script to run a differential accessibility analysis
## using DiffBind
## Sandra González de la Fuente
## Modified by María Santos Romero (28/10/2026)
##################################################################

## LOAD LIBRARIES ----

cat("\nChecking if all libraries are installed \n")
rm(list=ls())

# Cargar las librerías necesarias
suppressPackageStartupMessages({
library(DiffBind, quietly = TRUE)
library(BiocParallel, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(pheatmap, quietly = TRUE)
library(cluster, quietly = TRUE)
library(htmlwidgets, quietly = TRUE)
library(BiocManager, quietly = TRUE)
library(DESeq2, quietly = TRUE)
library(AnnotationDbi, quietly = TRUE)
library(scales)
library(ggrepel)
library(viridis)
library(colorspace)
library(Glimma)
library(optparse)
library(ChIPpeakAnno)
library(rtracklayer)
library(GenomicFeatures)
library(tidyverse)
})

## PARSE ARGUMENTS ----
option_list <- list(make_option(c("-d", "--design"), type="character", default=NULL, help="filename of design table", metavar="path"),
                    make_option(c("-l", "--log2FC"), type="numeric", default=0, help="log2FC", metavar="numeric"),
                    make_option(c("-g", "--gtf"), type="character", default=NULL, help="filename of gtf file", metavar="path"),
                    make_option(c("-o", "--output"), type="character", default=".", help="the name of the output directory", metavar="character"),
                    make_option(c("-c", "--cores"), type="integer", default=1, help="Number of cores", metavar="integer"),
                    make_option(c("-r", "--reference_condition"), type="character", default=NULL, help="Reference condition", metavar="character"))
                    
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$design)){
  print_help(opt_parser)
  stop("Please provide design table file.", call.=FALSE)
}
if (is.null(opt$log2FC)){
  print_help(opt_parser)
  stop("Please provide a log2FC.", call.=FALSE)
}
if (is.null(opt$reference_condition)){
  print_help(opt_parser)
  stop("Please provide a reference_condition", call.=FALSE)
}
if (is.null(opt$gtf)){
  print_help(opt_parser)
  stop("Please provide gtf file.", call.=FALSE)
}
if (is.null(opt$output)){
  print_help(opt_parser)
  stop("Please provide gtf file.", call.=FALSE)
}


## Assign input arguments to variables
reference_condition <- opt$reference_condition
design <- opt$design
log2fc_threshold <- opt$log2FC
out_dir <- opt$output
gtf_file <- opt$gtf

## !!FOR SCRIPT TESTING ----
#reference_condition <- "WT"
#design <- "/home/curso/dia3/diffbind.config"
#log2fc_threshold <- 1
#out_dir <-  "/home/curso/dia3/test_DiffBind"
#gtf_file <- "/home/curso/Software/dia3/reference/Mus_musculus.GRCm39.107.chr.gtf"

## PEAK TABLE ----
## create a csv file with SampleID, Condition, Replicate, bamReads Peaks Peakcaller PeakFormat, ScoreCol, Factor, Tissue
sampleDesign <- read.csv(design)

bamReads <- sampleDesign$bamReads
Peaks <- sampleDesign$Peaks
SampleID <- sampleDesign$SampleID

if(!all(names(bamReads) %in% SampleID)){
  stop("some names of bamReads is not in SampleID. names(bamReads)=",
       paste(names(bamReads), collapse = "; "),
       "SampleID=", paste(SampleID, collapse=";"))
}

names(Peaks) <- SampleID
rownames(sampleDesign) <- paste(sampleDesign$Condition, sampleDesign$Replicate, sep="_")
Condition <- sampleDesign$Condition
Replicate <- sampleDesign$Replicate
Peakcaller <- "macs3"
PeakFormat <- "narrowPeak"

samples <- data.frame(SampleID=SampleID,
                        Condition=Condition,
                        Replicate=Replicate,
                        bamReads=bamReads,
                        Peaks=Peaks,
                        Peakcaller=Peakcaller,
                        PeakFormat=PeakFormat,
                        ScoreCol=5)


if (!dir.exists(out_dir)){
  dir.create(out_dir)
}


# For each sample, read the peak calling files. Store it as a list; 
# list(sample1 = c(peak1, peak2, ...),
#      sample2 = c(peak1, peak2, peak3))
l <- lapply(Peaks, readLines)

# In case a peak calling file is empty, skip that sample
keep <- lengths(l)>0
samples <- samples[keep, , drop=FALSE]

# Sort samples based on condition and samples ID
samples <- samples[order(samples$Condition, samples$SampleID), , drop=FALSE]

if(nrow(samples)>3){
  
  write_csv(samples, file.path(out_dir, "samples.csv"))
  
  # Load your sample information and create an object for differential peak analysis
  db_obj <- dba(sampleSheet = file.path(out_dir, "samples.csv"))

  # Plot sample correlation based on called peaks
  png(file.path(out_dir, "DiffBind.sample.correlation.png"))
  plot(db_obj)
  dev.off()
  

  png(file.path(out_dir, "DiffBind.PCA.plot.png"))
  dba.plotPCA(db_obj,DBA_CONDITION, label=DBA_ID)
  dev.off()

  # Count sequencing reads in shared peaks regions
  db_obj <- dba.count(db_obj, bParallel=FALSE, summits=0, minOverlap=2, filter=5)

  # Visualizar el PCA basado en las condiciones y replicados
  pca_plot <- dba.plotPCA(db_obj, attributes=c(DBA_CONDITION), label=DBA_ID)
  
  png(file.path(out_dir, "DiffBind.PCA.plot.peakcount.png"))
  print(pca_plot)
  dev.off()
  
  png(file.path(out_dir, "DiffBind.sample.correlation.peakcount.png"))
  plot(db_obj)
  dev.off()
  
  # Save DiffBind object
  save(db_obj, file = file.path(out_dir, "DiffBind.RData"))
  gtf_file
  txdb <- txdbmaker::makeTxDbFromGFF(gtf_file)

  # Get ENSEMBL ID/symbol corresponde from the GTF
  gtf <- import(gtf_file)
  id2symbol <- function(gtf){
    if(is.null(gtf$gene_name)) return(NULL)
    x <- data.frame(id=gtf$gene_id, symbol=gtf$gene_name)
    x <- unique(x)
    x <- x[!duplicated(x$id), ]
    x <- x[!is.na(x$id), , drop=FALSE]
    if(nrow(x)==0) return(NULL)
    y <- x$symbol
    names(y) <- x$id
    y
  }

  id2symbol <- id2symbol(gtf)
  anno <- toGRanges(txdb)

  resList <- list()
  
  if(length(unique(Condition))>=2){

    contrasts <- combn(unique(Condition), m = 2, simplify = FALSE)
    names(contrasts) <- sapply(contrasts, function(.ele) paste(.ele[2], .ele[1], sep="-"))

    for(i in seq_along(contrasts)){
      gp1 <- as.character(contrasts[[i]][1])
      gp2 <- as.character(contrasts[[i]][2])

      if(sum(db_obj$samples$Condition %in% c(gp1, gp2))<3){
        next
      }

      db_obj <- dba.contrast(db_obj,
                            reorderMeta=list(Condition=reference_condition),
                            categories = DBA_CONDITION,
                            minMembers = 2)
    
      db_obj <- dba.normalize(db_obj, normalize="RLE", method=DBA_DESEQ2)
      
      pca_plot <- dba.plotPCA(db_obj, attributes=c(DBA_CONDITION), label=DBA_ID)
      png(file.path(out_dir, "DiffBind.PCA.plot.peakcount.norm.png"))
      print(pca_plot)
      dev.off()

      png(file.path(out_dir, "DiffBind.sample.correlation.peakcount.norm.png"))
      plot(db_obj)
      dev.off()
      
      # Perform statistical test for identifying regions with different accesibility
      db_obj <- dba.analyze(db_obj, bBlacklist = FALSE, bGreylist = FALSE,method=DBA_DESEQ2)

      # Retrieve results of differential accesibility analysis
      db_obj.DB <- dba.report(db_obj,method=DBA_DESEQ2, th=1, DataType=DBA_DATA_FRAME)
      db_obj.DB <- db_obj.DB[!is.na(db_obj.DB$Chr) &
                           !is.na(db_obj.DB$Start) &
                           !is.na(db_obj.DB$End), , drop=FALSE]

      db_obj.DB <- toGRanges(db_obj.DB)

      # Annotation
      db_obj.anno <- annotatePeakInBatch(db_obj.DB, AnnotationData = anno,
                                       output = "nearestLocation",
                                       PeakLocForDistance = "middle",
                                       FeatureLocForDistance = "TSS",
                                       ignore.strand = TRUE)

      if(length(id2symbol)>0) db_obj.anno$symbol[!is.na(db_obj.anno$feature)] <- id2symbol[db_obj.anno$feature[!is.na(db_obj.anno$feature)]]
      resList[[names(contrasts)[i]]] <- db_obj.anno[db_obj.anno$FDR<0.05]
      db_obj.m <- as.data.frame(unname(db_obj.anno),
                              stringsAsFactor=FALSE)

      # Export differential accesibility analysis results  
      write_tsv(db_obj.m, file.path(out_dir, paste0("DiffBind.res.", names(contrasts)[i], ".all.tsv")))
      
      # Filter results with FDR < 0.05
      db_obj.m <- db_obj.m[db_obj.m$FDR<0.05, ]
      write_tsv(db_obj.m, file.path(out_dir, paste0("DiffBind.res.", names(contrasts)[i], ".FDR.05.tsv")))
      
      # Filtrar los resultados para FDR < 0.05 y Fold Change absoluto > FC cutoof
      db_obj.m <- db_obj.m[db_obj.m$FDR < 0.05 & abs(db_obj.m$Fold) > log2fc_threshold, ]
      write_tsv(db_obj.m, file.path(out_dir, paste0("DiffBind.res.", names(contrasts)[i], ".FDR.05_log2FC.tsv")))

      ## VOLCANO PLOT ---
      # Prepare the data
      db_obj.m <- as.data.frame(unname(db_obj.anno),
                stringsAsFactor=FALSE)
      volcano_data <- as.data.frame(db_obj.m)
      volcano_data <- na.omit(volcano_data)  # Elimina filas con NA
      volcano_data <- volcano_data[!is.infinite(volcano_data$Fold) & !is.infinite(volcano_data$FDR), ]
      
      volcano_data$negLogFDR <- -log10(volcano_data$FDR)

      # Inicializa las columnas para las categorías y etiqueta
      volcano_data$Reg <- "NO"  # Inicializa con "NO"
      volcano_data$label <- NA   # Inicializa la columna de etiquetas como NA
      
      # Definir las condiciones para clasificar como UP o DOWN
      volcano_data$Reg[volcano_data$Fold > 0 & volcano_data$FDR < 0.05] <- "UP"    # Si Fold > 1 y FDR < 0.05
      volcano_data$Reg[volcano_data$Fold < -0 & volcano_data$FDR < 0.05] <- "DOWN"  # Si Fold < -1 y FDR < 0.05

      # Asignar las etiquetas a los puntos significativos
      volcano_data$label[volcano_data$Reg != "NO"] <-  volcano_data$symbol[volcano_data$Reg != "NO"]
      
      # Colores para los punto
      RegCol <- c("UP" = "#942025", "NO" = "grey", "DOWN" = "#788feb")
      regNames <- c("Downregulated", "Not significant", "Upregulated")

      # Define axis limits
      xlim_vals <- c(-4, 4)
      ylim_vals <- c(0, 10)

      # Define outliers as those points outside the axis limits
      volcano_data$Outlier <- (
        volcano_data$Fold < xlim_vals[1] |
        volcano_data$Fold > xlim_vals[2] |
        volcano_data$negLogFDR < ylim_vals[1] |
        volcano_data$negLogFDR > ylim_vals[2]
      )

      # Keep the original coordinates
      volcano_data$Fold_plot <- volcano_data$Fold
      volcano_data$FDR_plot <- volcano_data$negLogFDR

      # Place out-of-range points at the plot boundary
      volcano_data$Fold_plot <- pmax(xlim_vals[1], pmin(volcano_data$Fold_plot, xlim_vals[2]))
      volcano_data$FDR_plot <- pmax(ylim_vals[1], pmin(volcano_data$FDR_plot, ylim_vals[2]))
      
      volcano_plot <- ggplot(
        data = volcano_data,
        aes(x = Fold_plot, y = FDR_plot, col = Reg)
        ) +
        geom_point(
            aes(shape = Outlier),
            size = 2
        ) +
        theme_bw() +
        theme(
            panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
            panel.grid.minor = element_line(color = "grey90", linewidth = 0.3)
        ) +
        xlim(xlim_vals) +  # Límites del eje X
        ylim(ylim_vals) +  # Límites del eje Y
        geom_text_repel(
            aes(label = label),
            na.rm = TRUE,
            box.padding = 0.5,
            max.overlaps = 10
        ) +
        scale_color_manual(
            values = RegCol,
            labels = regNames,
            name = "Differentially\nAccessible Regions"
        ) +
        scale_shape_manual(
            values = c("FALSE" = 16, "TRUE" = 17),
            labels = c("FALSE" = "Regular", "TRUE" = "Outlier"),
            name = "Region type"
        ) +
        geom_hline(
            yintercept = -log10(0.05),
            linetype = "dashed",
            color = "black"
        )

        png(file.path(out_dir, paste0("DiffBind.", names(contrasts)[i], ".Volcano.plot.png")), width = 10, height = 10, units = "in", res = 300)
        print(volcano_plot)
        dev.off()

      ## MA PLOT ----
      # not normalized data
      png(file.path(out_dir, paste0("DiffBind.", names(contrasts)[i], ".MA.plot.png")))
      dba.plotMA(db_obj, bNormalize=F, method=DBA_DESEQ2)
      dev.off()

      # Normalized data
      png(file.path(out_dir, paste0("DiffBind.", names(contrasts)[i], "_norm.MA.plot.png")))
      dba.plotMA(db_obj, method=DBA_DESEQ2)
      dev.off()

      # Export counts table
      counts <- dba.peakset(db_obj, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
      write_tsv(counts, file.path(out_dir, paste0("DiffBind.", names(contrasts)[i], ".counts.tsv")))
      
      ##counts differential
      # Cambiar los nombres de las columnas de counts para que coincidan con db_obj.m
      diff_peaks <- db_obj.m[db_obj.m$FDR<0.05, ]
      colnames(counts)[2:3] <- c("start", "end")
    
      # Filtrar counts para incluir solo los picos diferenciales presentes en chip.
      differential_counts <- merge(diff_peaks <- db_obj.m[db_obj.m$FDR<0.05, ][, c("start", "end")], counts, by = c("start", "end"))

      ## PCA PLOT ----
      count_data <- differential_counts[, -c(1:3)]  # Eliminamos las columnas CHR, START, END

      # Transponer los datos para que las muestras sean las filas
      count_data_t <- t(count_data)

      # Realizar el PCA
      pca <- prcomp(count_data_t, scale. = TRUE)
      
      # Extraer las componentes principales
      pca_data <- as.data.frame(pca$x)
      
      # Añadir la columna 'Condition' a los resultados del PCA
      pca_data$Condition <- db_obj$samples$Condition
 
      
      ## HEATMAP_DIFF WITHOUT FC CUTOFF ----
      # Preparar los datos para el heatmap
      data_matrix <- differential_counts[, 4:ncol(differential_counts)]  # Seleccionar solo las columnas de datos de conteo
      data_matrix <- log2(data_matrix+1)

      # Escala de colores personalizada (granate para up, azul claro para down)
      hmcol <- colorRampPalette(c("lightblue", "white", "darkred"))(250)
      
      # Crear un data frame con la anotación de las condiciones
      df <- data.frame(Condition = db_obj$samples$Condition)
      rownames(df) <- db_obj$samples$SampleID
      colnames(data_matrix) <- rownames(df)
      
      # Definir colores para las condiciones
      conditions <- unique(db_obj$samples$Condition)
      subcolors <- list(Condition = setNames(c("darkgreen", "lightgreen"), conditions))

      #Calcular los tamaños dinámicos basados en el número de columnas y filas
      n <- pmax(ncol(data_matrix) - 6, 0)
      cw <- 20 - (5 / 21) * n  # Ajustar el tamaño de las celdas dinámicamente
      w <- 3.5 + (3.5 / 21) * n  # Ajustar el tamaño del gráfico
      if (nrow(data_matrix) < 2) {
        cr <- FALSE  # No hacer clustering en las filas si hay menos de 2 filas
      } else {
        cr <- TRUE}
      if (nrow(data_matrix) < 10) {
      h <- 4  # Ajustar la altura si hay menos de 10 filas
        } else {
      h <- 6  # Si hay más filas, usar una altura mayor
      }
      # Título para el gráfico
      titlePlot <- "Differential Peaks Heatmap"
      
      # Crear el heatmap con los parámetros especificados
      p <- pheatmap(
      mat = data_matrix,
      clustering_distance_rows="euclidean",
      cellwidth = cw,  # Tamaño dinámico de las celdas
      color = hmcol,  # Colores personalizados
      fontsize = 6.5,  # Tamaño de la fuente para los nombres de filas y columnas
      angle_col = 90,  # Ángulo de los nombres de las columnas
      fontsize_col = 7,  # Fuente para los nombres de las columnas
      border_color = FALSE,  # Sin borde para las celdas
      scale = "row",  # Escalar las filas
      cluster_cols = TRUE,  # Clustering de las columnas
      cluster_rows = cr,  # Clustering de las filas (depende del número de filas)
      annotation_col = df,  # Anotación de las columnas
      annotation_colors = subcolors,  # Colores para las condiciones
      show_rownames = FALSE,  # No mostrar nombres de las filas
      show_colnames = FALSE,  # Mostrar nombres de las columnas
      annotation_names_col = FALSE,  # No mostrar nombres de las columnas en la leyenda
      annotation_names_row = FALSE,  # No mostrar nombres de las filas en la leyenda
      main = paste0(titlePlot, "\n"),  # Título del gráfico
      treeheight_col = 10,  # Altura del árbol de las columnas
      annotation_legend = FALSE)  # No mostrar la leyenda de anotación)

      # Ajustar la visualización para que se vea bien en un formato PDF
      pdf(file.path(out_dir, paste0("DiffBind.", names(contrasts)[i], ".heatmapDiff.plot.pdf")),width = 10, height = 10,useDingbats = FALSE)
      print(p)
      dev.off()  # Cerrar el dispositivo gráfico

      ## HEATMAP_DIFF WITH FC CUTOFF ----
      # Cambiar los nombres de las columnas de counts para que coincidan con db_obj.m
      peakDiffLog <- db_obj.m[db_obj.m$FDR < 0.05 & abs(db_obj.m$Fold) > log2fc_threshold, ]
      colnames(counts)[2:3] <- c("start", "end")

      # Filtrar counts para incluir solo los picos diferenciales presentes en chip.
      differential_countslog <- merge(peakDiffLog[, c("start", "end")], counts, by = c("start", "end"))
      
      # Guardar el resultado si es necesario
           # Preparar los datos para el heatmap
      data_matrix <- differential_countslog[, 4:ncol(differential_countslog)]  # Seleccionar solo las columnas de datos de conteo
      
      # Escala de colores personalizada (granate para up, azul claro para down)
      data_matrix <- log2(data_matrix+1)
      hmcol <- colorRampPalette(c("lightblue", "white", "darkred"))(250)
      
      # Crear un data frame con la anotación de las condiciones
      df <- data.frame(Condition = db_obj$samples$Condition)
      rownames(df) <- db_obj$samples$SampleID
      colnames(data_matrix) <- rownames(df)
      
      # Extraer las condiciones de tu objeto de datos (por ejemplo, chip$samples$Condition)
      conditions <- unique(db_obj$samples$Condition)
      subcolors <- list(Condition = setNames(c("darkgreen", "lightgreen"), conditions))

      #Calcular los tamaños dinámicos basados en el número de columnas y filas
      n <- pmax(ncol(data_matrix) - 6, 0)
      cw <- 20 - (5 / 21) * n  # Ajustar el tamaño de las celdas dinámicamente
      w <- 3.5 + (3.5 / 21) * n  # Ajustar el tamaño del gráfico
      if (nrow(data_matrix) < 2) {
        cr <- FALSE  # No hacer clustering en las filas si hay menos de 2 filas
      } else {
        cr <- TRUE}
      if (nrow(data_matrix) < 10) {
        h <- 4  # Ajustar la altura si hay menos de 10 filas
      } else {
        h <- 6  # Si hay más filas, usar una altura mayor
      }
    
      # Crear el heatmap con los parámetros especificados
      p <- pheatmap(
        mat = data_matrix,
        clustering_distance_rows="euclidean",
        cellwidth = cw,  # Tamaño dinámico de las celdas
        color = hmcol,  # Colores personalizados
        fontsize = 6.5,  # Tamaño de la fuente para los nombres de filas y columnas
        angle_col = 90,  # Ángulo de los nombres de las columnas
        fontsize_col = 7,  # Fuente para los nombres de las columnas
        border_color = FALSE,  # Sin borde para las celdas
        scale = "row",  # Escalar las filas
        cluster_cols = TRUE,  # Clustering de las columnas
        cluster_rows = cr,  # Clustering de las filas (depende del número de filas)
        annotation_col = df,  # Anotación de las columnas
        annotation_colors = subcolors,  # Colores para las condiciones
        show_rownames = TRUE,  # No mostrar nombres de las filas
        show_colnames = TRUE,  # Mostrar nombres de las columnas
        annotation_names_col = FALSE,  # No mostrar nombres de las columnas en la leyenda
        annotation_names_row = FALSE,  # No mostrar nombres de las filas en la leyenda
        main = "Differential Peaks Heatmap\n",  # Título del gráfico
        treeheight_col = 10,  # Altura del árbol de las columnas
        annotation_legend = FALSE)  # No mostrar la leyenda de anotación)
      
      # Ajustar la visualización para que se vea bien en un formato PDF
      pdf(file.path(out_dir, paste0("DiffBind.", names(contrasts)[i], ".heatmapDiff_LFC.plot.pdf")),width = 10, height = 10,useDingbats = FALSE)
      print(p)
      dev.off()  # Cerrar el dispositivo gráfico
      
   
    }
  }

  resList <- if(length(resList)>1) GRangesList(resList) else if(length(resList)>0) resList[[1]]

  
  if(packageVersion("ChIPpeakAnno")>="3.23.12"){
    
    if(length(resList)>0){
      genomicElement <- genomicElementDistribution(resList,
                                        TxDb = txdb,
                                        promoterRegion=c(upstream=2000, downstream=500),
                                        geneDownstream=c(upstream=0, downstream=2000),
                                        promoterLevel=list(
                                          breaks = c(-2000, -1000, -500, 0, 500),
                                          labels = c("upstream 1-2Kb", "upstream 0.5-1Kb",
                                                     "upstream <500b", "TSS - 500b"),
                                          colors = c("#FFE5CC", "#FFCA99",
                                                     "#FFAD65", "#FF8E32")),
                                        plot = FALSE)
      genomicElement$plot <- genomicElement$plot +
        theme(
            plot.background = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA)
        )

      metagene_plot <- metagenePlot(resList, txdb)

      metagene_plot <- metagene_plot +
        labs(y = "Differential peaks") +
        theme(
            plot.background = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA)
        )

      #ggsave(file.path(out_dir, "genomicElementDistribuitonOfDiffBind.pdf"), plot=genomicElement$plot, width=9, height=9)
      ggsave(file.path(out_dir, "genomicElementDistribution_DiffBind.res.FDR.05.png"), plot=genomicElement$plot)
      
    #   ggsave(file.path(out_dir, "metagenePlotToTSSofDiffBind.pdf"), plot=out, width=9, height=9)
      ggsave(file.path(out_dir, "metagenePlotToTSS_DiffBind.res.FDR.05.png"), plot=metagene_plot)
    }

    if(length(samples$Peaks)>0){
      
      peaks <- mapply(samples$Peaks, samples$PeakFormat,
                      FUN=function(.ele, .format) toGRanges(.ele, format=.format),
                      SIMPLIFY = FALSE)

      names(peaks) <- samples$SampleID
      
      peaks <- GRangesList(peaks)
      
      genomicElement <- genomicElementDistribution(peaks,
                                        TxDb = txdb,
                                        promoterRegion=c(upstream=2000, downstream=500),
                                        geneDownstream=c(upstream=0, downstream=2000),
                                        promoterLevel=list(
                                          # from 5' -> 3', fixed precedence 3' -> 5'
                                          breaks = c(-2000, -1000, -500, 0, 500),
                                          labels = c("upstream 1-2Kb", "upstream 0.5-1Kb",
                                                     "upstream <500b", "TSS - 500b"),
                                          colors = c("#FFE5CC", "#FFCA99",
                                                     "#FFAD65", "#FF8E32")),
                                        plot = FALSE)
      genomicElement$plot <- genomicElement$plot +
        theme(
            plot.background = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA)
        )
      

      # metagene_plot <- metagenePlot(peaks, txdb)

      # metagene_plot <- metagene_plot +
      #   theme(
      #       plot.background = element_rect(fill = "white", colour = NA),
      #       panel.background = element_rect(fill = "white", colour = NA)
      #   )

      ggsave(file.path(out_dir, "genomicElementDistribution_DiffBind.allpeaks.png"), plot=genomicElement$plot)
      

     
      # ggsave(file.path(out_dir, "metagenePlotToTSS_DiffBind.allpeaks.png"), plot=metagene_plot)
      
      
      # Check peaks
      dup_peaks <- sum(unlist(lapply(peaks, duplicated)))
      cat("Duplicated ranges in peaks:", dup_peaks, "\n")
      
      # Check resList
      if(length(resList) > 0){
        if(class(resList) == "GRangesList"){
          dup_res <- sum(unlist(lapply(resList, duplicated)))
        } else {
          dup_res <- sum(duplicated(resList))
        }
        cat("Duplicated ranges in resList:", dup_res, "\n")
      }
      
      futile.logger::flog.threshold(
        futile.logger::ERROR,
        name = "VennDiagramLogger"
      )

      invisible(NULL)
      
      if(length(peaks)<=5){
        ol <- suppressMessages(findOverlapsOfPeaks(peaks))
        png(file.path(out_dir, "vennDiagram.all.png"))
        suppressMessages(makeVennDiagram(ol, connectedPeaks="keepAll"))
        dev.off()
   
      }else{
        peaks.split <- split(peaks, samples$Condition)
        peaks.split <- peaks.split[lengths(peaks.split)>1]
        null <- mapply(peaks.split, names(peaks.split), FUN=function(.peaks, .name){
          if(length(.peaks)<=5){
            ol <- invisible(findOverlapsOfPeaks(.peaks))
            png(file.path(out_dir, paste0("vennDiagram.", .name, ".png")))
            makeVennDiagram(ol, connectedPeaks="keepAll")
            dev.off()
          }
        })
      }
    }
  }
}

