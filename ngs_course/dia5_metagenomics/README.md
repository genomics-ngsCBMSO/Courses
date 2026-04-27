# DIA 5. Metagenomics - 16S

**Servicio de Análisis Biocomputacional (SABio); CBM** \
**Edition**: April, 2026 \
**Last update**: 27/04/2026

## Introduction to metagenomics

Metagenomics is a molecular biology and bioinformatics approach that enables the study of genetic material directly recovered from environmental samples, without the need to isolate or culture individual organisms. This is particularly important because the vast majority of microorganisms (estimated >99%) are not culturable using standard laboratory techniques.
By extracting and sequencing DNA directly from samples such as soil, water, or host-associated environments (e.g. gut, skin), metagenomics allows researchers to:

- Characterize microbial community composition
- Study microbial diversity and structure
- Infer ecological relationships
- Explore functional potential of microbial communities

### Amplicon-Based Metagenomics and 16S rRNA Sequencing

Amplicon metagenomics is a targeted sequencing approach, in contrast to shotgun metagenomics. Instead of sequencing all DNA present in a sample, a specific genetic marker is amplified using PCR and subsequently sequenced.
For bacteria and archaea, the most commonly used marker is the 16S ribosomal RNA (rRNA) gene, which:

- Is universally conserved across prokaryotes
- Contains both highly conserved and hypervariable regions (V1–V9)
- Allows taxonomic identification from phylum down to (often) genus level

In a 16S amplicon workflow:

Specific hypervariable regions of the 16S rRNA gene are amplified using primers
Amplicons are sequenced (typically Illumina paired-end)
Sequences are denoised and clustered into ASVs (Amplicon Sequence Variants)
ASVs are taxonomically annotated using reference databases
Diversity and community structure are analyzed

This repository implements a complete 16S rRNA amplicon analysis pipeline using **QIIME2**, a state-of-the-art platform for microbiome bioinformatics.

## Project Background

In this project, we analyze publicly available 16S rRNA gene sequencing data provided by the Schloss Lab as part of their well-established MiSeq Standard Operating Procedure (SOP) tutorial
(https://mothur.org/wiki/miseq_sop/).
The overarching biological question motivating this dataset is:

*How does natural variation in the gut microbiome relate to host physiology and stability over time?*

To address this question, the Schloss Lab conducted a longitudinal microbiome study in mice, focusing on the early post-weaning period, a critical developmental window characterized by rapid physiological and metabolic changes.

The study specifically investigates whether the rapid increase in body weight during the first ~10 days post-weaning is associated with changes in:

- Gut microbiome composition
- Community stability
- Temporal variability

These early dynamics are compared with the microbiome observed later in life (days 140–150), a period assumed to represent a more stable microbial ecosystem.

## 0. Environment setup

Load the conda environment "Dia1". This environment has installed a set of tools that we will use at different steps in the practice. For checking which software is available in this environment we can run `conda list`.


```bash
conda activate Dia1

# Create a new folder to store the data from today's practice
mkdir /home/curso/dia5
cd /home/curso/dia5
```

## 1. Data adquisition

### 1.1. Sequencing reads (FASTQ files)

Download the sequence reads that we’ll use in this analysis. These sequences are 250 bp and overlap in the V4 region of the 16S rRNA gene; this region is about 253 bp long. So looking at the files in the MiSeq_SOP folder that we’ve downloaded we will see 40 fastq files representing 10 time points from Female 3 and 1 mock community.


```bash
# Create a new folder to store fastq files
mkdir reads

#Illumina reads
wget https://www.mothur.org/w/images/d/d6/MiSeqSOPData.zip
unzip MiSeqSOPData.zip
mv MiSeq_SOP/*.fastq reads/
rm -rf MiSeq_SOP/
```

## 2.  Quality check with fastqc

Raw sequencing data (FASTQ files) is assessed for quality using tools like FastQC and MultiQC. This step detects issues such as adapter contamination, low-quality reads, or GC bias before continuing with further analysis.

```bash
mkdir fastqc
fastqc -o fastqc/ -t 6 reads/*.fastq
multiqc fastqc/* -o fastqc/
```

## 3.  Primer removal using cutadapt

Primer sequences used:

- Forward primer: CCTACGGGNGGCWGCAG
- Reverse primer: GACTACHVGGGTATCTAATCC

For each unique sample detected in the reads/ directory, trim the forward and reverse primers from its paired-end FASTQ files and save the cleaned reads.

```bash
mkdir cutadapt

for sample in `ls reads/*.fastq | cut -d"_" -f1 | sort -u`
do
    cutadapt -g CCTACGGGNGGCWGCAG -G GACTACHVGGGTATCTAATCC --cores 6 \
    -o cutadapt/${sample}_R1_trimmed.fastq \
    -p cutadapt/${sample}_R2_trimmed.fastq \
    reads/${sample}_*R1*.fastq reads/${sample}_*R2*.fastq
done
```

We use a loop to automatically process all samples. This works as follows:

1. `ls reads/*.fastq`  
Lists all FASTQ files in the reads/ directory.


2. `cut -d"_" -f1`   
Splits each filename using the underscore (_) as a delimiter and extracts the first field:

    ```bash
    F3D0_S188_L001_R1_001.fastq → F3D0  
    F3D0_S188_L001_R2_001.fastq → F3D0
    ```


3. `sort -u`  
Sorts the extracted sample names and removes duplicates, ensuring each sample is processed only once.


    The resulting list is used by the for loop, so the variable sample takes values such as:
    ```bash
    F3D0, F3D1, F3D2, ...
    ```

## 4. Manifest file generation

QIIME2 requires a manifest file mapping sample IDs to FASTQ file paths.

```bash
./make_manifest.sh -d cutadapt/
```
This produces manifest.csv, required for importing paired-end data.

## 5. Importing Data into Qiime2

All data that is used as input to QIIME 2 is in form of QIIME 2 artifacts, which contain information about the type of data and the source of the data. So, the first thing we need to do is import these sequence data files into a QIIME 2 artifact.

The semantic type of this QIIME 2 artifact is *SampleData[PairedEndSequencesWithQuality]*. *PairedEndSequencesWithQuality* QIIME 2 artifacts contain paired end sequences in fastq format (with quality) that are demultiplexed.

```bash
qiime tools import \
  --type SampleData[PairedEndSequencesWithQuality] \
  --input-path manifest.csv \
  --output-path paired_end_sequences.qza \
  --input-format PairedEndFastqManifestPhred33
```

## 6. Demultiplexing Summary

Always it’s useful to generate a summary of the reads imported. This allows us to determine how many sequences were obtained per sample, and also to get a summary of the distribution of sequence qualities at each position in our sequence data.

```bash
qiime demux summarize \
  --i-data paired_end_sequences.qza \
  --o-visualization paired_end_sequences.qzv
```

All QIIME 2 visualizers (i.e., commands that take a `--o-visualization parameter`) will generate a Visualization (i.e., a.qzv file). Visualizations can be viewed by loading them with QIIME 2 View or using the next command line:

```bash
qiime tools view paired_end_sequences.qzv
```

## 7. Denoising and ASV Inference with DADA2

**DADA2** is a pipeline for detecting and correcting (where possible) Illumina amplicon sequence data. As implemented in the q2-dada2 plugin, this quality control process will additionally filter any phiX reads (commonly present in marker gene Illumina sequence data) that are identified in the sequencing data, and will filter chimeric sequences.

The `dada2 denoise-paired` method requires two parameters that are used in quality filtering: `--p-trim-left m`, which trims off the first m bases of each sequence, and `--p-trunc-len n` which truncates each sequence at position n. This allows the user to remove low quality regions of the sequences. To determine what values to pass for these two parameters, we should review the Interactive Quality Plot tab in the `paired_end_sequences.qzv` file that was generated above.

```bash
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired_end_sequences.qza \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --o-table feature_table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats stats.qza \
  --verbose
```

This command generate QIIME 2 artifacts containing summary statistics. To view those summary statistics, we can visualize them using:

```bash
qiime metadata tabulate \
    --m-input-file stats.qza \
    --o-visualization stats.qzv

qiime tools view stats.qzv
```

## 8. FeatureTable and FeatureData summaries

After the quality filtering step completes, we’ll want to explore the resulting data. We can do this using the following two commands, which will create visual summaries of the data. The `feature-table summarize` command will give us information on how many sequences are associated with each sample and with each feature, histograms of those distributions, and some related summary statistics. The `feature-table tabulate-seqs` command will provide a mapping of feature IDs to sequences, and provide links to easily BLAST each sequence against the NCBI nt database.

We can add metadata information (biological, experimental, and contextual information)  using a **metadata file**. In our case, we will use the file *map_file.txt*.

Metadata play a central role in nearly all downstream analyses in QIIME2. Without a metadata file, sequencing data can only be analyzed in isolation.

QIIME2 enforces strict metadata validation rules to ensure reproducibility and correctness:

- The first column must contain sample IDs
- Sample IDs must exactly match those in the feature table
- Column names must be unique
- Missing values must be encoded as NA
- The file must be tab-separated, not comma-separated

```bash
qiime feature-table summarize \
  --i-table feature_table.qza \
  --o-visualization feature_table.qzv \
  --m-sample-metadata-file map_file.txt

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime tools view rep-seqs.qzv
```

## 9. Taxonomic Assignment Using BLAST

After denoising and inferring high‑resolution Amplicon Sequence Variants (ASVs), the next critical step in a 16S rRNA amplicon analysis is taxonomic assignment.
This step aims to answer the question:

**Which bacterial taxa do the inferred ASVs correspond to?**

Taxonomic assignment links exact nucleotide sequences (ASVs) to known microbial lineages, allowing biological interpretation of community composition.

In this pipeline, taxonomy is assigned using the *Greengenes 13_5* reference database. For that we need to import the reference sequences and taxonomy:

```bash
qiime tools import \
  --type FeatureData[Sequence] \
  --input-path gg_13_5.fasta \
  --output-path rep_seqs_gg.qza


qiime tools import \
  --type FeatureData[Taxonomy] \
  --input-path gg_13_5_taxonomy.txt \
  --output-path tax_gg.qza \
  --input-format HeaderlessTSVTaxonomyForm
```

The next step assigns taxonomy to each ASV using a **BLAST-based consensus classification** approach:

```bash
qiime feature-classifier classify-consensus-blast \
  --i-query rep-seqs.qza \
  --i-reference-taxonomy tax_gg.qza \
  --i-reference-reads rep_seqs_gg.qza \
  --o-search-results search_results.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
    --m-input-file search-results.qza \
    --o-visualization search_results.qzv

qiime metadata tabulate \
    --m-input-file taxonomy.qza \
    --o-visualization taxonomy.qzv

qiime tools view search_results.qzv
qiime tools view taxonomy.qzv
```

Next, we can view the taxonomic composition of our samples with interactive bar plots. Generate those plots with the following command and then open the visualization.

```bash
qiime taxa barplot \
  --i-table feature_table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file map_file.txt \
  --o-visualization taxa-bar-plots.qzv
```

We’re also often interested in performing a differential abundance test at a specific taxonomic level. To do this, we can collapse the features in our `FeatureTable[Frequency]` at the taxonomic level of interest. For example, we collapse our feature table at the genus level (i.e. level 6 of the Greengenes taxonomy).

```bash
qiime taxa collapse \
  --i-table feature_table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table collapsed_L6.qza

# This converts absolute read counts into relative abundances
qiime feature-table relative-frequency \
  --i-table collapsed_L6.qza \
  --o-relative-frequency-table collapsed_L6_relative.qza

# Converts a BIOM-format feature table into a human-readable tab-separated values (TSV) file
biom convert -i feature-table.biom -o feature-table.tsv --to-tsv
```

https://amplicon-docs.qiime2.org/en/stable/tutorials/moving-pictures/