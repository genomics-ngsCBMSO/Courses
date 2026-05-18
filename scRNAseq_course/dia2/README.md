# 🧬 Pipeline de Análisis Single-Cell RNA-seq (Día 2)

Este repositorio contiene los flujos de trabajo prácticos para la descarga/construcción de referencias genómicas, alineamiento y cuantificación con **Cell Ranger**, eliminación de ruido ambiental con **CellBender** e importación de datos en **R (Seurat)**.

## 🛠️ Paso 1: Gestión de Referencias Genómicas

### Opción A: Descarga de Referencia Oficial (Humano)
Ideal para organismos modelo estándar. Descarga directa desde los servidores de 10x Genomics.

```bash
# Crear y acceder al directorio de trabajo
mkdir -p reference && cd reference

# Descargar referencia GRCh38 oficial
wget "https://10xgenomics.com"

# Descomprimir el archivo tarball
tar -xzvf refdata-gex-GRCh38-2024-A.tar.gz
cd ..
```

### Opción B: Construcción de Referencia Personalizada (`cellranger mkref`)
Ejemplo práctico utilizando únicamente el Cromosoma Y obtenido desde Ensembl (Release 114).

```bash
mkdir -p chrY_example && cd chrY_example

# 1. Obtener y normalizar archivo FASTA (Secuencia)
wget https://ensembl.org
gzip -d Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz
mv Homo_sapiens.GRCh38.dna.chromosome.Y.fa chrY.fa

# Limpiar IDs de secuencia y pasar a mayúsculas para homogeneizar
awk '{if (\$1 ~ />/){print \(1} else {print\)_}}' chrY.fa | tr '[:lower:]' '[:upper:]' > chrY_mod.fa

# 2. Obtener y filtrar archivo GTF (Anotación)
wget https://ensembl.org
gzip -d Homo_sapiens.GRCh38.114.chr.gtf.gz

# Extraer metadatos iniciales y conservar únicamente las líneas asignadas al cromosoma Y
awk '\$1 ~ /^#/ {print \$0;next} {if (\$1 == "Y") print}' Homo_sapiens.GRCh38.114.chr.gtf > chrY.gtf

# 3. Indexar referencia con Cell Ranger
cellranger mkref --genome=chrY_index --fasta=chrY_mod.fa --genes=chrY.gtf
cd ..
```

---

## 🚀 Paso 2: Alineamiento y Cuantificación (`cellranger count`)

Mapeo de lecturas FASTQ, asignación de códigos de barras (Cell Barcodes), conteo de UMIs y generación de matrices de expresión digital.

```bash
mkdir -p cellranger && cd cellranger

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

cd ..
```

---

## 🧼 Paso 3: Corrección de ARN Ambiental (`cellbender`)

Aplicación de un modelo de aprendizaje profundo para remover background, contaminación por transcriptoma libre y falsas gotas positivas (*empty droplets*).

```bash
mkdir -p cellbender && cd cellbender

# Ejecutar remoción de ruido sobre la matriz cruda generada por Cell Ranger
cellbender remove-background \
           --input ../cellranger/PBMC/outs/raw_feature_bc_matrix.h5 \
           --output pbmc_cellbender_cleaned.h5

cd ..
```

---

## 📊 Paso 4: Carga y Análisis de Datos en R (`Seurat`)

Flujo final para importar las matrices procesadas dentro del ecosistema bioinformático de R e inicializar el objeto de análisis unicelular.

```R
# Iniciar sesión interactiva de R
R
```

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
