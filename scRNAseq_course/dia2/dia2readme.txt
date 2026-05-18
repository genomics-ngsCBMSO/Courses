## 1.2 Descarga referencias / Construir Referencias

mkdir reference
cd reference

#Podemos descargar la referencia ya hecha de la web oficial de cellranger (solo para humano, ratón y otros organismos modelo)

wget "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
tar -xzvf refdata-gex-GRCh38-2024-A.tar.gz

#O podemos construirla con cellranger mkref

mkdir chrY_example
cd chrY_example

wget https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz
gzip -d Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz
mv Homo_sapiens.GRCh38.dna.chromosome.Y.fa chrY.fa

#Quitamos todo lo que no sea ID de Secuencia
awk '{if ($1 ~ />/){print $1} else {print $_}}' chrY.fa | tr '[:lower:]' '[:upper:]' > chrY_mod.fa

wget https://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.chr.gtf.gz
gzip -d Homo_sapiens.GRCh38.114.chr.gtf.gz
awk '$1 ~ /^#/ {print $0;next} {if ($1 == "Y") print}' Homo_sapiens.GRCh38.114.chr.gtf > chrY.gtf


#Creamos la referencia con makeref

cellranger mkref --genome chrY --fasta chrY.fa --genes chrY.gtf 


########################
#DIA 2
########################

## 2.1 CELLRANGER: ALINEAMIENTO Y CUANTIFICACIÓN

mkdir cellranger
cd cellranger

ln -s ../dia1/reads/*gz ./

cellranger count --id=PBMC \
           --transcriptome=../dia1/reference/refdata-gex-GRCh38-2024-A/ \
           --fastqs=./ \
           --sample=pbmc_1k_v3 \
           --create-bam=true \
           --localcores=8 \
           --localmem=64

##2.2 CELLBENDER: FILTRADO DE RNA AMBIENTAL Y RUIDO

mkdir cellbender
cd cellbender

cellbender remove-background \
           --input raw_feature_bc_matrix.h5 \
           --output raw_feature_bc_matrix_cellbender.h5


## 2.2 LEER EN R

> R

library(Seurat)

data <- Read10X(data.dir = "filtered_feature_bc_matrix")
sp <- CreateSeuratObjetc(sp)
