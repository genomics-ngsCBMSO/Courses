#!/usr/bin/env Rscript
## R script to run a differential accessibility analysis
## using DiffBind
## Sandra González de la Fuente
##################################################################

## O. Libraries
############################################################################################################

cat("\n 0. Checking if all libraries are installed \n")
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
})


option_list <- list(make_option(c("-d", "--design"), type="character", default=NULL, help="filename of design table", metavar="path"),
                    make_option(c("-l", "--log2FC"), type="numeric", default=0, help="log2FC", metavar="numeric"),
                    make_option(c("-g", "--gtf"), type="character", default=NULL, help="filename of gtf file", metavar="path"),
                    make_option(c("-o", "--output"), type="character", default=".", help="the name of the output directory", metavar="character"),
                    make_option(c("-b", "--blacklist"), type="character", default=NULL, help="filename of blacklist file", metavar="path"),
                    make_option(c("-s", "--species"), type="character", default=NULL, help="species", metavar="string"),
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
out <- "sample.csv"
## create a csv file with SampleID, Condition, Replicate, bamReads Peaks Peakcaller PeakFormat, ScoreCol, Factor, Tissue
#register(MulticoreParam(workers = 24)) 
reference_condition <- opt$reference_condition
sampleDesign <- read.csv(opt$design)
#sampleDesign <- unique(sampleDesign[, !colnames(sampleDesign) %in% c("fastq_1", "fastq_2")])
log2fc_threshold <- opt$log2FC
bamReads <- sampleDesign$bamReads
#names(bamReads) <- sub(".mLb.clN.*bam", "", bamReads)
Peaks <- sampleDesign$Peaks
SampleID <- sampleDesign$SampleID
if(!all(names(bamReads) %in% SampleID)){
  stop("some names of bamReads is not in SampleID. names(bamReads)=",
       paste(names(bamReads), collapse = "; "),
       "SampleID=", paste(SampleID, collapse=";"))
}
#bamReads <- bamReads[SampleID]
names(Peaks) <- SampleID
rownames(sampleDesign) <- paste(sampleDesign$Condition, sampleDesign$Replicate, sep="_")
Condition <- sampleDesign$Condition
Replicate <- sampleDesign$Replicate
Peakcaller <- "macs2"
#PeakFormat <- sub("^.*?_peaks.(.*)$", "\\1", Peaks)
PeakFormat <- "narrowPeak"
block <- FALSE
if(any(grepl("treatment", colnames(sampleDesign), ignore.case = TRUE))){
  Treatment <- sampleDesign[SampleID, which(grepl("treatment",
                                                  colnames(sampleDesign),
                                                  ignore.case = TRUE))[1]]
  tt <- paste0(Condition, Treatment)
  if(length(unique(Treatment))>1 &&
     length(unique(tt))!=length(unique(Condition))){
    samples <- data.frame(SampleID=SampleID,
                          Condition=Condition,
                          Replicate=Replicate,
                          Treatment=Treatment,
                          bamReads=bamReads,
                          Peaks=Peaks,
                          Peakcaller=Peakcaller,
                          PeakFormat=PeakFormat,
                          ScoreCol=5)
    block <- TRUE
  }else{
    samples <- data.frame(SampleID=SampleID,
                          Condition=Condition,
                          Replicate=Replicate,
                          bamReads=bamReads,
                          Peaks=Peaks,
                          Peakcaller=Peakcaller,
                          PeakFormat=PeakFormat,
                          ScoreCol=5)
  }
}else{
  samples <- data.frame(SampleID=SampleID,
                        Condition=Condition,
                        Replicate=Replicate,
                        bamReads=bamReads,
                        Peaks=Peaks,
                        Peakcaller=Peakcaller,
                        PeakFormat=PeakFormat,
                        ScoreCol=5)
}

pf <- opt$output
dir.create(pf)
l <- lapply(Peaks, readLines)
keep <- lengths(l)>0
samples <- samples[keep, , drop=FALSE]
samples <- samples[order(samples$Condition, samples$SampleID), , drop=FALSE]

BLACKLIST <- paste0("DBA_BLACKLIST_", toupper(opt$species))
if(exists(BLACKLIST)){
  BLACKLIST <- get(BLACKLIST)
}else{
  BLACKLIST <- FALSE
}

if(nrow(samples)>3){
  write.csv(samples, file.path(pf, "sample.csv"))

  chip <- dba(sampleSheet = file.path(pf, "sample.csv"))
  png(file.path(pf, "DiffBind.sample.correlation.png"))
  plot(chip)
  dev.off()
  png(file.path(pf, "DiffBind.PCA.plot.png"))
  dba.plotPCA(chip,DBA_CONDITION, label=DBA_ID)
  dev.off()

  chip <- dba.count(chip, bParallel=FALSE,summits=0,minOverlap=2,filter=5)
  save(chip, file = file.path(pf, "chip.RData"))
  #load("chip_machos.RData")

  chip <- dba.blacklist(chip, blacklist=BLACKLIST, greylist=FALSE)
  chip.bk <- chip

  # Visualizar el PCA basado en las condiciones y replicados
  pca_plot <- dba.plotPCA(chip, attributes=c(DBA_CONDITION), label=DBA_ID)
  png(file.path(pf, "PCA_plot_samples_count.png"))
  print(pca_plot)
  dev.off()
  png(file.path(pf, "PCA_plot_samples_count_correlation.png"))
  plot(chip)
  dev.off()
  txdb <- txdbmaker::makeTxDbFromGFF(opt$gtf)
  gtf <- import(opt$gtf)
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
      if(sum(chip.bk$samples$Condition %in% c(gp1, gp2))<3){
        next
      }
      chip <- chip.bk
      if(!block){
        chip <- dba.contrast(chip,
                             reorderMeta=list(Condition=reference_condition),
                             categories = DBA_CONDITION,minMembers = 2)
      }else{
        chip <- dba.contrast(chip,
                             reorderMeta=list(Condition=reference_condition),
                             categories = DBA_CONDITION,
                             block = DBA_TREATMENT,minMembers = 2)
      }
      #chip <- dba.normalize(chip,background=TRUE, library="DBA_LIBSIZE_BACKGROUND",method=DBA_DESEQ2)
      chip <- dba.normalize(chip, normalize="RLE",method=DBA_DESEQ2)
      #save(chip, file = file.path(pf, "chipNorm.RData"))
      #load("chipNorm.RData")
      pca_plot <- dba.plotPCA(chip, attributes=c(DBA_CONDITION), label=DBA_ID)
      png(file.path(pf, "PCA_plot_samples_norm.png"))
      print(pca_plot)
      dev.off()
      png(file.path(pf, "PCA_plot_samples_norm_correlation.png"))
      plot(chip)
      dev.off()
      chip <- dba.analyze(chip, bBlacklist = FALSE, bGreylist = FALSE,method=DBA_DESEQ2)
      chip.DB <- dba.report(chip,method=DBA_DESEQ2, th=1, DataType=DBA_DATA_FRAME)
      chip.DB <- chip.DB[!is.na(chip.DB$Chr) &
                           !is.na(chip.DB$Start) &
                           !is.na(chip.DB$End), , drop=FALSE]
      chip.DB <- toGRanges(chip.DB)

      # Annotation
      chip.anno <- annotatePeakInBatch(chip.DB, AnnotationData = anno,
                                       output = "nearestLocation",
                                       PeakLocForDistance = "middle",
                                       FeatureLocForDistance = "TSS",
                                       ignore.strand = TRUE)
      if(length(id2symbol)>0) chip.anno$symbol[!is.na(chip.anno$feature)] <- id2symbol[chip.anno$feature[!is.na(chip.anno$feature)]]
      resList[[names(contrasts)[i]]] <- chip.anno[chip.anno$FDR<0.05]
      chip.m <- as.data.frame(unname(chip.anno),
                              stringsAsFactor=FALSE)
      write.csv(chip.m, file.path(pf, paste0("DiffBind.res.", names(contrasts)[i], ".all.csv")))
      chip.m <- chip.m[chip.m$FDR<0.05, ]
      write.csv(chip.m, file.path(pf, paste0("DiffBind.res.", names(contrasts)[i], ".FDR.05.csv")))
      # Filtrar los resultados para FDR < 0.05 y Fold Change absoluto > 0.5
      chip.m <- chip.m[chip.m$FDR < 0.05 & abs(chip.m$Fold) > log2fc_threshold, ]
      write.csv(chip.m, file.path(pf, paste0("DiffBind.res.", names(contrasts)[i], ".FDR.05_log2FC.csv")))
      # plots
      pdf(file.path(pf, paste0("DiffBind.", names(contrasts)[i], ".PCA-contrast.plot.pdf")))
      dba.plotPCA(chip, contrast=1, method=DBA_DESEQ2,attributes=DBA_CONDITION, label=DBA_ID)
      dev.off()
      
      ##volcano
      chip.m <- as.data.frame(unname(chip.anno),
                stringsAsFactor=FALSE)
      volcano_data <- as.data.frame(chip.m)
      volcano_data <- na.omit(volcano_data)  # Elimina filas con NA
      volcano_data <- volcano_data[!is.infinite(volcano_data$Fold) & !is.infinite(volcano_data$FDR), ]
      #rownames(volcano_data) <- volcano_data$symbol
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
      # Ajustar los límites para los ejes X e Y
      #xlim_vals <- c(-4, 4)
      #ylim_vals <- c(0, 25)
      py <- -log10(volcano_data$FDR)
      ylim_vals <- c(min(py),max(py))
      px <- volcano_data$Fold
      xlim_vals <- c(min(px),max(px)) 
      volcano_plot <- ggplot(data = volcano_data, aes(x = Fold, y = -log10(FDR), col = Reg, label = label)) +
      geom_point() +  # Dibujar los puntos
      theme_minimal() +  # Tema limpio sin cuadrícula
      theme(panel.grid = element_blank()) +  # Eliminar la cuadrícula
      geom_text_repel(aes(label = label), na.rm = TRUE,box.padding = 0.5, max.overlaps = 10) +  # Etiquetas para puntos significativos
      scale_color_manual(values = RegCol, labels = regNames, name = "Differentially\nAccesible Regions") +  # Colores y etiquetas de la leyenda
      xlim(xlim_vals) +  # Límites del eje X
      ylim(ylim_vals) +  # Límites del eje Y
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")  # Línea horizontal de corte en p-value = 0.05
      pdf(file.path(pf, paste0("DiffBind.", names(contrasts)[i], ".Volcano_Diff.plot.pdf")), width = 10, height = 10, useDingbats = FALSE)
      print(volcano_plot)
      dev.off()
      png(file.path(pf, paste0("DiffBind.", names(contrasts)[i], ".MA.plot.png")))
      dba.plotMA(chip, method=DBA_DESEQ2)
      dev.off()
      png(file.path(pf, paste0("DiffBind.", names(contrasts)[i], ".MA.plot_notNorm.png")))
      dba.plotMA(chip, bNormalize=F,method=DBA_DESEQ2)
      dev.off()

      png(file.path(pf, paste0("DiffBind.", names(contrasts)[i], ".Volcano.plot.png")))
      tryCatch(dba.plotVolcano(chip,  method=DBA_DESEQ2,bUsePval = FALSE), error=function(.e) message(.e))
      dev.off()

      # export counts table
      counts <- dba.peakset(chip, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
      write.csv(counts, file.path(pf, paste0("DiffBind.", names(contrasts)[i], ".counts.csv")))
      
      ##counts differential
      # Cambiar los nombres de las columnas de counts para que coincidan con chip.m
      chipDiff <- chip.m[chip.m$FDR<0.05, ]
      colnames(counts)[2:3] <- c("start", "end")
      # Filtrar counts para incluir solo los picos diferenciales presentes en chip.
      differential_counts <- merge(chipDiff[, c("start", "end")], counts, by = c("start", "end"))
      # Guardar el resultado si es necesario

      ##pca DIFF
      # Datos simulados
      
      count_data <- differential_counts[, -c(1:3)]  # Eliminamos las columnas CHR, START, END
      # Transponer los datos para que las muestras sean las filas
      count_data_t <- t(count_data)
      # Realizar el PCA
      pca <- prcomp(count_data_t, scale. = TRUE)
      # Extraer las componentes principales
      pca_data <- as.data.frame(pca$x)
      # Añadir la columna 'Condition' a los resultados del PCA
      pca_data$Condition <- chip$samples$Condition
      # Crear el gráfico PCA con colore
      pca_plot <- ggplot(pca_data, aes(PC1, PC2, fill = Condition, color = Condition)) +
      geom_point(size = 4, shape = 21) +  # Puntos en el gráfico
      geom_text(aes(label = rownames(pca_data)), size = 3, hjust = -0.1, vjust = -0.5) +  # Etiquetas de muestras
      ggtitle("Principal Components Plot") +
      theme_light() +
      scale_y_continuous(
      paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 2), "%)"),
      expand = c(0.05, 0)) +
      scale_x_continuous(
      paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 2), "%)"),
      expand = c(0.05, 0)) +
      scale_fill_manual(values = c("darkgreen", "lightgreen")) +  # Colores para las condiciones
      scale_color_manual(values = c("darkolivegreen", "forestgreen")) + # Colores para las condiciones
      theme(
      plot.title = element_text(hjust = 0.5, size = 10, margin = margin(0, 0, 10, 0)),
      axis.title.x = element_text(size = 9, margin = margin(10, 0, 0, 0)),
      axis.title.y = element_text(size = 9, margin = margin(0, 10, 0, 0)),
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.15, "in"))      
      pdf(file.path(pf, paste0("DiffBind.", names(contrasts)[i], ".PCA_Diff.plot.pdf")),width = 7, height = 5,useDingbats = FALSE)
      print(pca_plot)     
      dev.off()
      
      ##HEATMAP_DIFF
      # Preparar los datos para el heatmap
      data_matrix <- differential_counts[, 4:ncol(differential_counts)]  # Seleccionar solo las columnas de datos de conteo
      data_matrix <- log2(data_matrix+1)
      # Escala de colores personalizada (granate para up, azul claro para down)
      hmcol <- colorRampPalette(c("lightblue", "white", "darkred"))(250)
      # Crear un data frame con la anotación de las condiciones
      df <- data.frame(Condition = chip$samples$Condition)
      rownames(df) <- chip$samples$SampleID
      colnames(data_matrix) <- rownames(df)
      # Definir colores para las condiciones
      conditions <- unique(chip$samples$Condition)
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
      show_colnames = TRUE,  # Mostrar nombres de las columnas
      annotation_names_col = FALSE,  # No mostrar nombres de las columnas en la leyenda
      annotation_names_row = FALSE,  # No mostrar nombres de las filas en la leyenda
      main = paste0(titlePlot, "\n"),  # Título del gráfico
      treeheight_col = 10,  # Altura del árbol de las columnas
      annotation_legend = FALSE)  # No mostrar la leyenda de anotación)
      # Ajustar la visualización para que se vea bien en un formato PDF
      pdf(file.path(pf, paste0("DiffBind.", names(contrasts)[i], ".heatmapDiff.plot.pdf")),width = 10, height = 10,useDingbats = FALSE)
      print(p)
      dev.off()  # Cerrar el dispositivo gráfico

      ##counts differential log2FC
      # Cambiar los nombres de las columnas de counts para que coincidan con chip.m
      chipDiffLog <- chip.m[chip.m$FDR < 0.05 & abs(chip.m$Fold) > log2fc_threshold, ]
      colnames(counts)[2:3] <- c("start", "end")
      # Filtrar counts para incluir solo los picos diferenciales presentes en chip.
      differential_countslog <- merge(chipDiffLog[, c("start", "end")], counts, by = c("start", "end"))
      # Guardar el resultado si es necesario
           # Preparar los datos para el heatmap
      data_matrix <- differential_countslog[, 4:ncol(differential_countslog)]  # Seleccionar solo las columnas de datos de conteo
      # Escala de colores personalizada (granate para up, azul claro para down)
      data_matrix <- log2(data_matrix+1)
      hmcol <- colorRampPalette(c("lightblue", "white", "darkred"))(250)
      # Crear un data frame con la anotación de las condiciones
      df <- data.frame(Condition = chip$samples$Condition)
      rownames(df) <- chip$samples$SampleID
      colnames(data_matrix) <- rownames(df)
      # Extraer las condiciones de tu objeto de datos (por ejemplo, chip$samples$Condition)
      conditions <- unique(chip$samples$Condition)
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
      show_rownames = TRUE,  # No mostrar nombres de las filas
      show_colnames = TRUE,  # Mostrar nombres de las columnas
      annotation_names_col = FALSE,  # No mostrar nombres de las columnas en la leyenda
      annotation_names_row = FALSE,  # No mostrar nombres de las filas en la leyenda
      main = paste0(titlePlot, "\n"),  # Título del gráfico
      treeheight_col = 10,  # Altura del árbol de las columnas
      annotation_legend = FALSE)  # No mostrar la leyenda de anotación)
      # Ajustar la visualización para que se vea bien en un formato PDF
      pdf(file.path(pf, paste0("DiffBind.", names(contrasts)[i], ".heatmapDiff_LFC.plot.pdf")),width = 10, height = 10,useDingbats = FALSE)
      print(p)
      dev.off()  # Cerrar el dispositivo gráfico
      
      #log2FC
      # Extraer resultados de DiffBind
      chip.m <- as.data.frame(unname(chip.anno),
                stringsAsFactor=FALSE)
      volcano_data <- as.data.frame(chip.m)
      volcano_data <- na.omit(volcano_data)  # Elimina filas con NA
      volcano_data <- volcano_data[!is.infinite(volcano_data$Fold) & !is.infinite(volcano_data$FDR), ]
      volcano_data$Reg <- "NO"  # Inicializa con "NO"
      volcano_data$label <- NA   # Inicializa la columna de etiquetas como NA
      # Definir las condiciones para clasificar como UP o DOWN
      volcano_data$Reg[volcano_data$Fold > log2fc_threshold & volcano_data$FDR < 0.05] <- "UP"    # Si Fold > 1 y FDR < 0.05
      volcano_data$Reg[volcano_data$Fold < -log2fc_threshold & volcano_data$FDR < 0.05] <- "DOWN"  # Si Fold < -1 y FDR < 0.05
      # Asignar las etiquetas a los puntos significativos
      volcano_data$label[volcano_data$Reg != "NO"] <-  volcano_data$symbol[volcano_data$Reg != "NO"]
      # Colores para los punto
      RegCol <- c("UP" = "#942025", "NO" = "grey", "DOWN" = "#788feb")
      regNames <- c("Downregulated", "Not significant", "Upregulated")
      py <- -log10(volcano_data$FDR)
      ylim_vals <- c(min(py),max(py))
      px <- volcano_data$Fold
      xlim_vals <- c(min(px),max(px)) 
      volcano_plot <- ggplot(data = volcano_data, aes(x = Fold, y = -log10(FDR), col = Reg, label = label)) +
      geom_point() +  # Dibujar los puntos
      theme_minimal() +  # Tema limpio sin cuadrícula
      theme(panel.grid = element_blank()) +  # Eliminar la cuadrícula
      geom_text_repel(aes(label = label), na.rm = TRUE,box.padding = 0.5, max.overlaps = 10) +  # Etiquetas para puntos significativos
      scale_color_manual(values = RegCol, labels = regNames, name = "Differentially\nAccesible Regions") +  # Colores y etiquetas de la leyenda
      xlim(xlim_vals) +  # Límites del eje X
      ylim(ylim_vals) +  # Límites del eje Y
      geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), colour = "black", linewidth=0.2, linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")  # Línea horizontal de corte en p-value = 0.05
      
      
      pdf(file.path(pf, paste0("DiffBind.Volcano_Diff_LFC", names(contrasts)[i], ".pdf")),width = 10, height = 10,useDingbats = FALSE)
      print(volcano_plot)
      dev.off()
   
    }
  }
  resList <- if(length(resList)>1) GRangesList(resList) else if(length(resList)>0) resList[[1]]

  if(packageVersion("ChIPpeakAnno")>="3.23.12"){
    if(length(resList)>0){
      out <- genomicElementDistribution(resList,
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

      ggsave(file.path(pf, "genomicElementDistribuitonOfDiffBind.pdf"), plot=out$plot, width=9, height=9)
      ggsave(file.path(pf, "genomicElementDistribuitonOfDiffBind.png"), plot=out$plot)
      out <- metagenePlot(resList, txdb)
      ggsave(file.path(pf, "metagenePlotToTSSofDiffBind.pdf"), plot=out, width=9, height=9)
      ggsave(file.path(pf, "metagenePlotToTSSofDiffBind.png"), plot=out)
    }
    if(length(samples$Peaks)>0){
      peaks <- mapply(samples$Peaks, samples$PeakFormat,
                      FUN=function(.ele, .format) toGRanges(.ele, format=.format),
                      SIMPLIFY = FALSE)
      names(peaks) <- samples$SampleID
      peaks <- GRangesList(peaks)
      out <- genomicElementDistribution(peaks,
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

      ggsave(file.path(pf, "genomicElementDistribuitonOfEachPeakList.pdf"), plot=out$plot, width=9, height=9)
      ggsave(file.path(pf, "genomicElementDistribuitonOfEachPeakList.png"), plot=out$plot)

      out <- metagenePlot(peaks, txdb)
      ggsave(file.path(pf, "metagenePlotToTSSOfEachPeakList.pdf"), plot=out, width=9, height=9)
      ggsave(file.path(pf, "metagenePlotToTSSOfEachPeakList.png"), plot=out)
      
      
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
      
      
      if(length(peaks)<=5){
        ol <- findOverlapsOfPeaks(peaks)
        png(file.path(pf, "vennDiagram.all.png"))
        makeVennDiagram(ol, connectedPeaks="keepAll")
        dev.off()
      }else{
        peaks.split <- split(peaks, samples$Condition)
        peaks.split <- peaks.split[lengths(peaks.split)>1]
        null <- mapply(peaks.split, names(peaks.split), FUN=function(.peaks, .name){
          if(length(.peaks)<=5){
            ol <- findOverlapsOfPeaks(.peaks)
            png(file.path(pf, paste0("vennDiagram.", .name, ".png")))
            makeVennDiagram(ol, connectedPeaks="keepAll")
            dev.off()
          }
        })
      }
    }
  }
}

