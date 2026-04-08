#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(BiocManager, quietly = TRUE)
  library(clusterProfiler, quietly = TRUE)
  library(enrichplot, quietly = TRUE)
  library(UpSetR, quietly = TRUE)
  library(ggplot2, quietly = TRUE)
  library(RColorBrewer, quietly = TRUE)
  library(optparse, quietly = TRUE)
  library(cowplot, quietly = TRUE)
  library(ggridges, quietly = TRUE)
  library(UpSetR, quietly = TRUE)
  
##Barplot specific

library(stringr)
library(dplyr)
library(forcats)
})

# Parámetros para el script
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Archivo de entrada con genes y valores [ID, stat]", metavar="character"),
  make_option(c("-o", "--orgdb"), type="character", default="org.Mm.eg.db", help="OrgDb para GSEA [default: %default]", metavar="character"),
  make_option(c("-p", "--prefix"), type="character", default="fosfoproteinas_16E_vs_cntrl_all", help="Prefijo para los archivos de salida [default: %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Debe proveer un archivo de entrada con la opción -i", call.=FALSE)
}

orgdb <- opt$orgdb
prefix <- opt$prefix
input <- opt$input

# Carga de la base de datos de genes
suppressPackageStartupMessages({
  library(orgdb, quietly = TRUE, character.only = TRUE)
})

# Lectura de datos
data <- read.table(input, sep= "\t", header=T, quote="")
data <- data[, c("ID", "log2FC", "qvalue")]

#data$stat <- data$log2FC*-log10(data$pvalue)

data <- data[, c("ID", "log2FC")]
# Configuración de GSEA
key <- "ENSEMBL"
ont <- "ALL"
terms <- NULL

if (!is.null(terms)){
  terms <- unlist(strsplit(terms, ",")) 
} else {
  terms <- NULL
}

# Preparación de los datos para GSEA
colnames(data) <- c("ID", "log2FC")
dat <- data$log2FC
names(dat) <- as.character(data$ID)
dat_nofcNA <- dat[!is.na(dat)]
genes <- names(dat_nofcNA)
dat_sort <- sort(dat_nofcNA, decreasing=TRUE)


# Realización de GSEA con diferentes valores de p-value y una semilla fija para reproducibilidad
set.seed(1234)
#for (pvalue in c(0.05, 0.1, 0.5)) {
  #egs <- gseGO(geneList = dat_sort, OrgDb = orgdb, keyType = key, ont = ont, verbose = T, pvalueCutoff = pvalue, seed = T, eps = 0)
  #print(paste("egs Returned with:", length(egs@result$ID), "IDS at cutoff pvalue", pvalue))
#}

# Definir el p-value final
pvalue <- 0.5
egs <- gseGO(geneList = dat_sort, OrgDb = orgdb, keyType = key, ont = ont, verbose = T, pvalueCutoff = pvalue, seed = T, eps = 0)

# Convertir los resultados a nombres legibles
egs_genename <- setReadable(egs, OrgDb = orgdb)

# Guardar la sesión
save.image(file = "sessionSave.RData")

# Guardar resultados en archivos
write.table(egs, file = paste("tableGO_", prefix, ".txt", sep =""), sep= "\t", quote = F)
write.table(egs_genename, file = paste("tableGO_", prefix, "_genename.txt", sep =""), sep= "\t", quote = F)

# Graficar Dotplot
tiff(file = paste("GO_", prefix, "_dotplot.tiff", sep =""), units = 'in', width = 10, height = 10, res = 200)
dotplot(egs, x = "GeneRatio", color = "p.adjust", showCategory = 15, font.size = 10, split = ".sign", label_format = 50) + facet_grid(.~.sign)
dev.off()

tiff(file = paste("GO_", prefix, "_dotplot_pvalue.tiff", sep =""), units = 'in', width = 10, height = 10, res = 200)
dotplot(egs, x = "GeneRatio", color = "pvalue", showCategory = 15, font.size = 10, split = ".sign", label_format = 50) + facet_grid(.~.sign)
dev.off()

# Simplificar resultados
result_simplify <- simplify(egs)
gene_count <- egs@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
egs_filt <- left_join(egs@result, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)
egs_filt$sign <- ifelse(egs_filt$NES > 0, "activated", "suppressed")
egs_filt_subset <- head(egs_filt, 50)

# Graficar Barplot
tiff(file = paste("GO_", prefix, "_barplots.tiff", sep =""), units = 'in', width = 10, height = 10, res = 200)
ggplot(egs_filt_subset, aes(x = NES, y = fct_reorder(Description, NES))) +
  geom_col(aes(fill = GeneRatio), color = "black", width = 0.6) +
  scale_fill_gradient(limits = c(0, 1), low = "green", high = "red") +
  labs(x = "Normalized Enrichment Score (NES)", y = "Gene Ontology Term") +
  theme_classic()
dev.off()

# Simplificación de los resultados
gene_count <- result_simplify@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
egs_filt <- left_join(result_simplify@result, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)
egs_filt$sign <- ifelse(egs_filt$NES > 0, "activated", "suppressed")
egs_filt_subset <- head(egs_filt, 20)

# Graficar Barplot simplificado
tiff(file = paste("GO_", prefix, "_barplots_simplify.tiff", sep =""), units = 'in', width = 10, height = 10, res = 200)
ggplot(egs_filt_subset, aes(x = NES, y = fct_reorder(Description, NES))) +
  geom_col(aes(fill = GeneRatio), color = "black", width = 0.6) +
  scale_fill_gradient(limits = c(0, 1), low = "green", high = "red") +
  labs(x = "Normalized Enrichment Score (NES)", y = "Gene Ontology Term") +
  theme_classic()
dev.off()

# Modificar color del Barplot
pdf(file = paste("GO_", prefix, "_barplots_selected_colour_v2.pdf", sep = ""), useDingbats = FALSE, width = 10, height = 10)
egs_filt_subset$NES_direction <- ifelse(egs_filt_subset$NES < 0, "Negative", "Positive")
ggplot(egs_filt_subset, aes(x = NES, y = fct_reorder(Description, NES))) +
  geom_col(aes(fill = ifelse(NES < 0, -GeneRatio, GeneRatio)), color = "black", width = 0.6) +
  scale_fill_gradientn(colors = c("#244ee3", "#d2daf7", "white", "white", "#fcd4d5", "#942025"), values = c(0, 0.25, 0.49, 0.51, 0.75, 1), name = "GeneRatio", breaks = c(-1, -0.5, 0, 0.5, 1), na.value = "white") +
  labs(x = "Normalized Enrichment Score (NES)", y = "Gene Ontology Term") +
  theme_classic()
dev.off()

# Gene-concept network
tiff(file = paste("GO_", prefix, "_gene_concept_net.tiff", sep =""), units = 'in', width = 10, height = 10, res = 200)
cnetplot(egs_genename, categorySize = "NES", foldChange = dat_sort, font.size = 10, showCategory = 5)
dev.off()

# Enrichment map
a <- pairwise_termsim(egs, method = "JC", semData = NULL, showCategory = 20)
tiff(file = paste("GO_", prefix, "_enrich_map.tiff", sep =""), units = 'in', width = 10, height = 10, res = 200)
emapplot(a, showCategory = 15, color = "NES")
dev.off()

# Ridgeline plot
tiff(file = paste("GO_", prefix, "_ridge.tiff", sep =""), units = 'in', width = 15, height = 15, res = 200)
ridgeplot(egs) + labs(x = "Expression fold change") + theme_classic()
dev.off()

##Upset plot (of the 20 first terms)

genes_20first <- as.data.frame(as.factor(head(egs@result$core_enrichment, 20)))
lista_20first <- list()
for (i in 1:nrow(genes_20first)){
    lista_20first[[i]] <- unlist(strsplit(as.character(genes_20first[i,1]),split="/"))   
}
uniq_genes <- as.character(unique(names(dat_sort)))

func_20first <- egs$Description[1:20]
mat <- matrix(0L, nrow = length(uniq_genes), ncol = length(func_20first)) 

for (i in 1:length(uniq_genes)) {
for (j in 1:length(func_20first)) {
gen <- uniq_genes[i]
if (gen %in% lista_20first[[j]]) {
mat[i,j] =  1
}}} 
 
mat_20first <- as.data.frame(mat)
colnames(mat_20first) <- func_20first
row.names(mat_20first) <- uniq_genes


jpeg(file = paste(prefix, "_upset_20first.jpeg", sep =""), units = 'in', width = 15, height = 10, res = 300)
    upset(mat_20first, nsets=10, order.by="freq", sets.bar.color="skyblue")
invisible(dev.off())
