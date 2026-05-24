# DIA 1. scRNA-seq DATA DOWNLOAD AND READ QUALITY ASSESMENT
**Servicio de Análisis Biocomputacional (SABio); CBM** \
**Edition**: May, 2026 \
**Last update**: 22/05/2026

## Setup

Create a folder for `Day1`practice and load the conda environment with tools for scRNAseq data analyisis installed.

```bash
# Create a folder
mkdir /home/curso/dia1

# Load the conda environment
conda activate sc-RNAseq
```

## Download scRNAseq reads

Reads were downloaded from 10X Genomics Datasets (https://www.10xgenomics.com/datasets) a repository with public or example  sequencing datasets generated using 10x Genomics technologies, especially:

- Single-cell RNA sequencing (scRNA-seq)
- Single-cell ATAC-seq (chromatin accessibility)
- Spatial transcriptomics
- Multiome (RNA + chromatin from same cell)

The repository includes data generated for different tissues and organisms. This data is commonly used for training or pipeline/tools development or setup.

The data for this practice corresponds to peripheral blood mononuclear cells (PBMCs) from a healthy human donor. 

**Data source**: https://www.10xgenomics.com/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0

Downloaded files: 

| File name | Content | Sample |
|---        |---      | --- |
|pbmc_1k_v3_S1_L001_I1_001.fastq.gz | Sample barcode | Sample01 | 
pbmc_1k_v3_S1_L002_I1_001.fastq.gz | Sample barcode | Sample01 |
pbmc_1k_v3_S1_L001_R1_001.fastq.gz | Cell + UMI barcodes | Sample01 |
pbmc_1k_v3_S1_L002_R1_001.fastq.gz | Cell + UMI barcodes | Sample01 |
pbmc_1k_v3_S1_L001_R2_001.fastq.gz | 3' RNA | Sample01 |
pbmc_1k_v3_S1_L002_R2_001.fastq.gz | 3' RNA | Sample01 |

*Note*: for the course the data has already been downloaded and can be found in /home/curso/Software/curso_scRNAseq


```bash
# Create a directory to store the reads
mkdir /home/curso/dia1/reads
cd /home/curso/dia1/reads

# Create a symbolic link to the read FATSQ files
ln -s /home/curso/Software/curso_scRNAseq/reads/*.fastq.gz .
```

Explore the content and structure of the read files

```bash
# Show the first lines  
    # from an index file
    zcat pbmc_1k_v3_S1_L001_I1_001.fastq.gz | head 

    # from an R1 file
    zcat pbmc_1k_v3_S1_L001_R1_001.fastq.gz | head 

    # from an R2 file
    zcat pbmc_1k_v3_S1_L001_R2_001.fastq.gz | head 
   
# Calculate the number of lines per FASTQ
for file in *.fastq.gz
do 
    echo ${file}
    zcat ${file} | wc -l
done

# Calculate the number of reads per file
for file in *.fastq.gz
do 
    echo ${file}
    lines=`zcat ${file} | wc -l`
    echo "Numer of reads: $((lines / 4))" # $((...)) evaluates what is inside as a mathematical expression 
done

# Get the read length for each file
for file in *.fastq.gz
do 
    echo ${file}
    read_length=`zcat ${file} | awk 'NR % 4 == 2 { print length($0) }' | head -n 1`
    echo ${read_length}
done

# Print the 10 most prevalent sequences in the index file
for file in *_I1_001.fastq.gz
do 
    zcat ${file} | awk 'NR%4==2' | sort | uniq -c | sort -nr | head
done
```

## Assess read quality

Run `FASTQC` to assess reads quality

```bash
# Create a dir to store quality reports
mkdir /home/curso/dia1/fastqc

# Run FASTQC
fastqc /home/curso/dia1/reads/*.fastq.gz -t 6 --outdir /home/curso/dia1/fastqc

# Aggregate reports using multiqc
cd /home/curso/dia1/fastqc
/home/curso/miniconda3/bin/multiqc .
```