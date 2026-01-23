#!/usr/bin/env Rscript


### ENVIROMENT AND PARAMETERS SETUP
# Load libraries
suppressPackageStartupMessages({
  library(BiocManager, quietly = TRUE)
  library(clusterProfiler, quietly = TRUE)
  library(enrichplot, quietly = TRUE)
  library(ggplot2, quietly = TRUE)
  library(RColorBrewer, quietly = TRUE)
  library(optparse, quietly = TRUE)
  library(cowplot, quietly = TRUE)
  library(ggridges, quietly = TRUE)
  library(stringr)
  library(dplyr)
  library(forcats)
})

# Define input parameters
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Archivo de entrada con genes y valores. Debe contener las columnas 'ID' y 'stat' [ID, stat]", metavar="character"),
  make_option(c("-o", "--orgdb"), type="character", default="org.Mm.eg.db", help="OrgDb para GSEA [default: %default]", metavar="character"),
  make_option(c("-p", "--prefix"), type="character", default="fosfoproteinas_16E_vs_cntrl_all", help="Prefijo para los archivos de salida [default: %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Stop if input file was not provided
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Debe proveer un archivo de entrada con la opción -i", call.=FALSE)
}

# Parse parameters
orgdb <- opt$orgdb
prefix <- opt$prefix
input <- opt$input

# Load specific ORGB database
suppressPackageStartupMessages({
  library(orgdb, quietly = TRUE, character.only = TRUE)
})


### DATA LOADING
# Read the input data table
input_table <- read.table(input, sep= "\t", header=T, quote="")

# Select ID and stat collumns
data <- input_table[, c("ID", "stat")]

# Format input data to be used with GSEA
stat_values <- data$stat
names(stat_values) <- as.character(data$ID)

# Remove entries with no stat value
stat_values_noNA <- stat_values[!is.na(stat_values)]

# Get the list of gene IDs after filtering
genes <- names(stat_values_noNA)

# Sort the data based on the stat value
data_sort <- sort(stat_values_noNA, decreasing=TRUE)


### RUN GSEA
# GSEA configuration parameters
key <- "ENSEMBL" # input ID key
ont <- "ALL" # terms ontology
pvalue <- 0.05 # p-value threshold for enriched terms
terms <- NULL # custom list of terms (optional)

# Parse custom terms, if any
if (!is.null(terms)){
  terms <- unlist(strsplit(terms, ",")) 
} else {
  terms <- NULL
}

# Run GSEA
set.seed(1234) # seed is used to ensure reproducibility 
egs <- gseGO(geneList = data_sort, OrgDb = orgdb, keyType = key, ont = ont, verbose = T, pvalueCutoff = pvalue, seed = T, eps = 0)

# Change gene IDs to gene symbols
egs_genename <- setReadable(egs, OrgDb = orgdb)

# Simplify results to reduce the number of redundant terms
result_simplify <- simplify(egs)

# Count genes in each enriched GO term
gene_count <- egs@result %>% 
  group_by(ID) %>%  # group rows by GO term ID
  summarise(count = sum(str_count(core_enrichment, "/")) + 1) # column containing genes contributing to enrichment (e.g BRCA/mTOR/...)

# Add gene count an compute gene ratio
egs_annotated <- left_join(egs@result, gene_count, by = "ID") %>% 
              mutate(GeneRatio = count/setSize) # setSize = total number of genes in the GO term

# Classify the enrichment direction
egs_annotated$sign <- ifelse(egs_annotated$NES > 0, "activated", "suppressed")

# Keep the top 20 GO terms
egs_top <- head(egs_annotated, 20)


### EXPORT RESULTS
# Save R session and GSEA results object
save.image(file = "sessionSave.RData")
saveRDS(egs, file = "gsea_results.rds")

# Save results as a table
# with ENSEMBL IDs
write.table(egs, file = paste("tableGO_", prefix, ".txt", sep =""), sep= "\t", quote = F)

# with gene symbols
write.table(egs_genename, file = paste("tableGO_", prefix, "_genename.txt", sep =""), sep= "\t", quote = F)


### GENERATE RESULTS PLOTS
# Define a list to store the plots
gsea_plot <- list()

# Generate result plots
# dotplot with padj
gsea_plot$dotplot_padj <- dotplot(egs, x = "GeneRatio", color = "p.adjust", showCategory = 15, font.size = 10, split = ".sign", label_format = 50) + facet_grid(.~.sign)

# dotplot with pval
gsea_plot$dotplot_pval <- dotplot(egs, x = "GeneRatio", color = "pvalue", showCategory = 15, font.size = 10, split = ".sign", label_format = 50) + facet_grid(.~.sign) 

# barplot
gsea_plot$barplot <- ggplot(egs_top, aes(x = NES, y = fct_reorder(Description, NES))) +
                      geom_col(aes(fill = GeneRatio), color = "black", width = 0.6) +
                      scale_fill_gradient(limits = c(0, 1), low = "green", high = "red") +
                      labs(x = "Normalized Enrichment Score (NES)", y = "Gene Ontology Term") +
                      theme_classic()

# barplot top terms
gsea_plot$barplot_top <- ggplot(egs_top, aes(x = NES, y = fct_reorder(Description, NES))) +
                          geom_col(aes(fill = GeneRatio), color = "black", width = 0.6) +
                          scale_fill_gradient(limits = c(0, 1), low = "green", high = "red") +
                          labs(x = "Normalized Enrichment Score (NES)", y = "Gene Ontology Term") +
                          theme_classic()

# barplot enrichment direction
gsea_plot$barplot_direction <- ggplot(egs_top, aes(x = NES, y = fct_reorder(Description, NES))) +
                                geom_col(aes(fill = ifelse(NES < 0, -GeneRatio, GeneRatio)), color = "black", width = 0.6) +
                                scale_fill_gradientn(colors = c("#244ee3", "#d2daf7", "white", "white", "#fcd4d5", "#942025"), values = c(0, 0.25, 0.49, 0.51, 0.75, 1), name = "GeneRatio", breaks = c(-1, -0.5, 0, 0.5, 1), na.value = "white") +
                                labs(x = "Normalized Enrichment Score (NES)", y = "Gene Ontology Term") +
                                theme_classic()

# Gene concept network
gsea_plot$cnetplot <- cnetplot(egs_genename, showCategory = 20, foldChange=data_sort)

# String plot
sim_matrix <- pairwise_termsim(egs)
gsea_plot$emapplot <- emapplot(sim_matrix, showCategory = 20, color = "pvalue")

# Ridgeline plot
gsea_plot$ridgeplot <- ridgeplot(egs) + labs(x = "Expression fold change") + theme_classic()

# Save all plot
for (plot_name in names(gsea_plot)){
    plot <- gsea_plot[[plot_name]]
    # as png
    ggsave(plot, file=paste0(plot_name,"_", ont, "_", prefix, ".png"), device = "png", bg="white", height = 10, width = 10)
    # as pdf
    ggsave(plot, file=paste0(plot_name,"_", ont, "_", prefix, ".pdf"), device = "pdf", height = 10, width = 10)
}
