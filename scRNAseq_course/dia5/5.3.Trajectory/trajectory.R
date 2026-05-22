#SCRIPT PARA REALIZAR ANALISIS DE TRAYECTORIAS A PARTIR DE SCRNASEQ

###############################
#0. CARGAR LIBRERIAS Y SEMILLA
###############################

library(Seurat)
library(openxlsx)
library(dplyr)
library(patchwork)
library(SingleCellExperiment)
library(monocle3)
library(ggplot2)

set.seed(1234)

#Cargamos el nuevo dataset con datos adecuados para analisis de trayectoria
cds <- readRDS("trayectory_dataset_CDS.rds")

cds

#Podemos ver las metadatos tanto de las filas (genes) como de las columnas (células)
rowData(cds)
colData(cds)

#########################################
#2. PROCESAMIENTO
#########################################

#Reducimos dimensionalidad (UMAP)
cds <- reduce_dimension(cds)

jpeg("UMAP.jpeg", units = "in", height = 10, width = 12, res = 300)
plot_cells(cds, 
           cell_size = 0.5) +
           theme(axis.text = element_text(size = 15),
                 axis.title = element_text(size = 15))
dev.off()

#Como la anotación ya está hecha podemos colorear las células por tipo celular

plot_cells(cds, 
           cell_size = 0.5,  
           color_cells_by = "cell.type",
           group_label_size = 5) +
           theme(axis.text = element_text(size = 15),
                 axis.title = element_text(size = 15))

#Monocle recalcula clusters internamente
cds <- cluster_cells(cds)

jpeg("Clusters.jpeg", units = "in", height = 10, width = 12, res = 300)
plot_cells(cds,
           cell_size = 0.5,
           color_cells_by = "cluster",
           group_label_size = 7) +
           theme(axis.text = element_text(size = 15),
                 axis.title = element_text(size = 15))
dev.off()

#########################################
#3. ANALISIS DE TRAYECTORIA
#########################################

#3.1 Construye un grafo que conecta células (estructura de trayectoria)
cds <- learn_graph(cds)

jpeg("Trayectory_celltype.jpeg", units = "in", height = 10, width = 12, res = 300)
plot_cells(cds,
           color_cells_by = "cell.type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           group_label_size= 3)
dev.off()     

#3.2 Definir el nodo raíz

#El pseudotiempo necesita un punto inicial biológico
#Ejemplo: elegir manualmente un cluster o tipo celular inicial

#Visualiza primero para decidir
plot_cells(
  cds,
  color_cells_by = "cell.type",
  label_groups_by_cluster = TRUE,
  graph_label_size=5
)

#Ejemplo: seleccionar células de tipo "Monocytes" como raíz
#root_cells <- colnames(cds)[colData(cds)$celltype == "Monocytes"]
#cds <- order_cells(cds, root_cells = root_cells)

#Alternativa: Seleccionar nodo raiz manualmente
cds <- order_cells(cds)


#3.3 VISUALIZACION DE PSEUDOTIEMPO

jpeg("Pseudotime.jpeg", units = "in", height = 10, width = 12, res = 300)
plot_cells(
  cds,
  cell_size = 0.5,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)
dev.off()

#########################################
#4. GENES QUE CAMBIAN EN PSEUDOTIEMPO
#########################################

#Detectar genes dinámicos a lo largo de la trayectoria (Tarda un poco)
deg_pseudotime <- graph_test(cds, neighbor_graph = "principal_graph")

#Ver los más significativos
head(deg_pseudotime)

#Visualizar genes dinámicos
#Seleccionar genes top
top_genes <- rownames(subset(deg_pseudotime, q_value < 0.05))[1:3

plot_cells(
  cds,
  genes = top_genes,
  show_trajectory_graph = TRUE
)

#########################################
#12. GUARDAR RESULTADOS
#########################################

saveRDS(cds, "monocle3_trajectory.rds")
