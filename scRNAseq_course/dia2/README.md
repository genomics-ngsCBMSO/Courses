# 🧬 Pipeline de análisis Single-cell RNA-seq (Día 2)

Este repositorio contiene los flujos de trabajo prácticos para la descarga/construcción de referencias genómicas, alineamiento y cuantificación con **Cell Ranger**, eliminación de ruido ambiental con **CellBender** e importación de datos en **R (Seurat)**.

## 🛠️ Paso 1: Gestión de referencias genómicas

### Opción A: Descarga de Referencia oficial (Humano)
Ideal para organismos modelo estándar. Descarga directa desde los servidores de 10x Genomics.

```bash
# Crear y acceder al directorio de trabajo
mkdir reference
cd reference

# Descargar referencia GRCh38 oficial
wget "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"

# Descomprimir el archivo tarball
tar -xzvf refdata-gex-GRCh38-2024-A.tar.gz
```

### Opción B: Construcción de Referencia personalizada (cellranger mkref)
Ejemplo práctico utilizando únicamente el cromosoma Y obtenido desde Ensembl (Release 114).

```bash
mkdir chrY_example
cd chrY_example

# 1. Obtener y normalizar archivo FASTA (Secuencia)
wget https://ensembl.org
gzip -d Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz
mv Homo_sapiens.GRCh38.dna.chromosome.Y.fa chrY.fa

# Limpiar IDs de secuencia y pasar a mayúsculas para homogeneizar
awk '{if (\$1 ~ />/){print \(1} else {print\)_}}' chrY.fa | tr '[:lower:]' '[:upper:]' > chrY_mod.fa

# 2. Obtener y filtrar archivo GTF (Anotación)
wget https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz
gzip -d Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz
mv Homo_sapiens.GRCh38.dna.chromosome.Y.fa chrY.fa
awk '{if ($1 ~ />/){print $1} else {print $_}}' chrY.fa | tr '[:lower:]' '[:upper:]' > chrY_mod.fa #Quitamos todo lo que no sea ID de Secuencia

wget https://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.chr.gtf.gz
gzip -d Homo_sapiens.GRCh38.114.chr.gtf.gz
awk '$1 ~ /^#/ {print $0;next} {if ($1 == "Y") print}' Homo_sapiens.GRCh38.114.chr.gtf > chrY.gtf #Quitamos todo lo que no sea del chrY

# 3. Indexar referencia con Cell Ranger
cellranger mkref --genome chrY --fasta chrY.fa --genes chrY.gtf
```

---

## 🚀 Paso 2: Alineamiento y Cuantificación (`cellranger count`)

Mapeo de lecturas FASTQ, asignación de códigos de barras (Cell Barcodes), conteo de UMIs y generación de matrices de expresión digital.

```bash
mkdir cellranger
cd cellranger

# Crear enlaces simbólicos a las lecturas del Día 1 para optimizar espacio en disco
ln -s ../dia1/reads/*gz ./

# Ejecutar el pipeline de conteo cuantitativo
cellranger count --id=PBMC \
           --transcriptome=../reference/refdata-gex-GRCh38-2024-A/ \
           --fastqs=./ \
           --sample=pbmc_1k_v3 \
           --create-bam=true \
           --localcores=8 \
           --localmem=64
```

---

## 🧼 Paso 3: Corrección de ARN ambiental (cellbender)

Aplicación de un modelo de aprendizaje profundo para eliminar background, contaminación por transcriptoma libre y falsas gotas positivas (*empty droplets*).

```bash
mkdir cellbender
cd cellbender

# Ejecutar cellbender
cellbender remove-background \
           --input ../cellranger/PBMC/outs/raw_feature_bc_matrix.h5 \
           --output pbmc_cellbender_cleaned.h5
```

---

## 📊 Paso 4: Carga y análisis de datos en R (Seurat)

Flujo final para importar las matrices procesadas dentro de R e inicializar el objeto de análisis.


```R
# Cargar librería analítica
library(Seurat)

# 1. Leer los directorios MEX filtrados generados
data_dir <- "cellranger/PBMC/outs/filtered_feature_bc_matrix"
counts <- Read10X(data.dir = data_dir)

# 2. Inicializar el objeto Seurat básico
pbmc_object <- CreateSeuratObject(counts = counts, project = "PBMC_1K", min.cells = 3, min.features = 200)

# Visualizar estructura del objeto creado
print(pbmc_object)
```
