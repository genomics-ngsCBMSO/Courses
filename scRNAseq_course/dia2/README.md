# 🧬 Single-cell RNA-seq Analysis steps (Day 2)

This repository contains the practical workflows for downloading/constructing genomic references, alignment and quantification with **Cell Ranger**, ambient RNA removal with **CellBender**, and data import into **R (Seurat)**.

## 🛠️ Step 1: Genomic reference management

### Option A: Official Reference download (Human)
Ideal for standard model organisms. Direct download from 10x Genomics servers.

```bash
# Crear y acceder al directorio de trabajo
mkdir reference
cd reference

# Descargar referencia GRCh38 oficial
wget "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"

# Descomprimir el archivo tarball
tar -xzvf refdata-gex-GRCh38-2024-A.tar.gz
```

### Option B: Custom Reference construction (`cellranger mkref`)
Practical example using only chromosome Y obtained from Ensembl (Release 114).

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
awk '{if (\$1 ~ />/){print \(1} else {print\)_}}' chrY.fa | tr '[:lower:]' '[:upper:]' > chrY_mod.fa #Quitamos todo lo que no sea ID de Secuencia

wget https://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.chr.gtf.gz
gzip -d Homo_sapiens.GRCh38.114.chr.gtf.gz
awk '\$1 ~ /^#/ {print \$0;next} {if (\$1 == "Y") print}' Homo_sapiens.GRCh38.114.chr.gtf > chrY.gtf #Quitamos todo lo que no sea del chrY

# 3. Indexar referencia con Cell Ranger
cellranger mkref --genome chrY --fasta chrY.fa --genes chrY.gtf
```

---

## 🚀 Step 2: Alignment and Quantification (`cellranger count`)

Mapping of FASTQ reads, cell barcode assignment, UMI counting, and digital expression matrix generation.

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

## 🧼 Step 3: Ambient RNA correction (`cellbender`)

Application of a deep learning model to eliminate background, free-transcriptome contamination, and false-positive droplets (*empty droplets*).

```bash
mkdir cellbender
cd cellbender

# Ejecutar cellbender
cellbender remove-background \
           --input ../cellranger/PBMC/outs/raw_feature_bc_matrix.h5 \
           --output raw_feature_bc_matrix_cellbender.h5
```

As the previous result follows the initial training trend but suffers from large fluctuations and a drastic drop at the end, and furthermore does not predict the cell number well, we will apply a smaller --learning-rate and add the command: --expected-cells

```bash
# Ejecutar cellbender cambiando parámetros
cellbender remove-background \
cellbender remove-background --input raw_feature_bc_matrix.h5   --expected-cells 1221 y --learning-rate 0.000005 
--output raw_feature_bc_matrix_cellbender_ECLR.h5

```

---

## 📊 Step 4: Data loading and analysis in R (`Seurat`)

Final workflow to import the processed matrices into R and initialize the analysis object.

```R
# Cargar librería analítica
library(Seurat)

# 1. Leer los directorios MEX filtrados generados
data_dir <- "cellranger/PBMC/outs/filtered_feature_bc_matrix"
data <- Read10X(data.dir = data_dir)

# 2. Inicializar el objeto Seurat básico
seu <- CreateSeuratObject(counts = data)

# Visualizar estructura del objeto creado
head(seu)
```
