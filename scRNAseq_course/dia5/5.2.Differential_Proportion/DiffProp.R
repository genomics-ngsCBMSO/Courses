#SCRIPT PARA REALIZAR ANALISIS DE PROPORCIONES DE TIPOS CELULAREs

###############################
#0. CARGAR LIBRERIAS Y SEMILLA
###############################

library(Seurat)
library(openxlsx)
library(dplyr)
library(ggplot2)

set.seed(1234)

#En este caso necesitamos comparar entre muestras.
#Se ha creado un objeto con la combinación de nuestros datos (Donor1) y datos de otra muestra (Donor2)

seu <- readRDS("combined_donors.rds")

seu
head(seu)

#Mirar cuantas celulas tenemos para cada donor

table(seu$Sample)

#################################
#1. Extraer los datos del objeto
#################################

#Columna con cell type -> predicted.id
#Columna con donante -> Sample
    
md <- seu@meta.data %>%
      select(CellType = predicted.id, Sample = Sample)

#Las dos columnas deben ser factores
md$CellType <- factor(md$CellType)
md$Sample <- factor(md$Sample)

head(md)

##########################################
#2. Frecuencias y test estadístico global
##########################################

#Calculamos las frecuencias de cada tipo celular por muestra
freqs <- table(md$Sample, md$CellType)
freqs

#Test estadístico de frecuencias entre los dos donantes: X cuadrado

#Dado que en datos de scRNA-seq es común tener tipos celulares poco abundantes (frecuencias bajas), 
#usamos simulate.p.value = TRUE para calcular el p-value mediante simulaciones en lugar de la aproximación teórica, 
#lo que produce resultados más robustos en tablas de contingencia dispersas

chisq_res <- chisq.test(freqs, simulate.p.value=TRUE)
chisq_res

#p-value = 0.0004998
#Hay diferencias significativas en las proporciones de tipos celulares entre los dos donantes

##########################################
#3. Test estadístico individual
##########################################

#Cargamos la función ya hecha para calcular test estadistico individual

source("ctprop_test.R")

#Calculamos el porcentaje de cada tipo celular en cada donante

freqs_df <- as.data.frame(freqs)
colnames(freqs_df) <- c("Sample", "CellType", "Frequency")

freqs_df$Proportion <- freqs_df$Frequency / ave(freqs_df$Frequency, freqs_df$Sample, FUN = sum)

#Aplicamos la función y guardamos los resultados

results <- ctprop_test(freqs_df, freqs)

write.xlsx(results, "DiffProp_results.xlsx", row.names = FALSE)

##############################
#4. Plots de los resultados
##############################

#Barplot

jpeg("Barplot_Prop.jpeg", units = "in", height = 10, width = 12, res = 300)
ggplot(freqs_df, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = scales::percent(Proportion, accuracy = 0.1)),
            position = position_stack(vjust = 0.5),
            size = 3, fontface = "bold") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Proporción de tipos celulares por muestra",
       y = "Proporción",
       x = "Sample") +
  theme_minimal(base_size = 15)
dev.off()
  
#PieChart

jpeg("PieChar_Prop.jpeg", units = "in", height = 10, width = 12, res = 300)
ggplot(freqs_df, aes(x = "", y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", width = 1) +
  geom_text(aes(label = scales::percent(Proportion, accuracy = 0.1)),
            position = position_stack(vjust = 0.5),
            size = 3, fontface = "bold") +
  coord_polar(theta = "y") +
  facet_wrap(~ Sample) +
  theme_void(base_size = 15) +
  labs(title = "Composición celular por muestra")
dev.off()



  
