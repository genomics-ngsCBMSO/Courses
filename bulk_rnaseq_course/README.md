# INTRODUCTION

This mini-course on bulk RNA-seq analysis aims to teach the basis of the bulk RNA-seq workflow and to develop the skills required to run analyses on a High-Performance Computing (HPC) system.

The training covers basic Linux commands, Conda environments, HPC usage, and RNA-seq analysis using publicly available data from a zebrafish experiment comparing wild-type and mutant samples. 

## Basic Linux commands

```bash
ls                                  # list directory contents
ls -ltrh                            # long format, sorted by time, reverse order, human-readable sizes

mkdir dir_name                      # create one directory
mkdir dir1_name dir2_name           # create multiple directories
mkdir -p project/results            # create nested directories

pwd                                 # showing your current path 

cd                                  # change directory
cd dir_name                         # go into a directory
cd ..                               # go up one level
cd ~                                # go to home directory
cd -                                # go back to previous directory

rm                                  # remove files (!permanent)
rm file.txt                         # remove file
rm -d dir_name                      # remove empty directory
rm -r dir_name                      # remove directory recursively

mv file.txt new_destination         # move file
mv oldname.txt  newname.txt         # rename file
mv *.txt folder/                    # move multiple files

wget                                # download files from the web
wget https://example.com/file.zip   # download file

gzip file.txt                       # compress file to file.txt.gz
gzip -d file.txt.gz                 # decompress file

echo "Hello"                        # print text
echo $PATH                          # print a variable
tree                                # show full directory tree

cat file.txt                        # print entire file
less file.txt                       # scrollable view (best for large files)
head file.txt                       # first 10 lines
tail file.txt                       # last 10 lines
```

## File formats
    
**FASTA**:

The first line starts with ">" followed by the description or identifier of the sequence.
The following line is the actual sequence itself in standard one-letter code.

*Example*:

    >Sequence_id
    NUCLEOTIDESEQUENCE

    >chr1
    AAATATTATGCGCGAGTTTCAGAAA


**FASTQ**:

Represents nucleotide sequneces, but also contains the correspondig quality of each nucleotide. Four lines per sequence:
- Line 1 begins with '@' followed by sequence identifier and an optional description.
- Line 2 is the raw sequence of nucleotides.
- Line 3 begins with a '+' and is a free text field.
- Line 4 encodes the quality values for the sequence in Line 2. Quality (phred score) = probability that the corresponding basecall is incorrect.

*Example*:

    @SEQ_ID
    GAGAGTTTGCGAGCTTTGCTAGCT
    +
    !''*-(((***+))%%%++)-%%%

**SAM**:

Text-based format for storing sequences aligned to a reference sequence. \
Consists of a header and an alignment section. \
In brief, it consist of a header section and reads (with other information) in tab delimited format.


**BAM**:

Binary representation of a SAM file. \
Compressed so no human-readable.
    
**GFF/GFF3/GTF**:

General feature format (gff) is used for describing genes and other features od DNA, RNA and protien sequences.
Standard annotation of genomes.

*Example*:
    
    ##description: evidence-based annotation of the human genome (GRCh38), version 25 (Ensembl 85)
    ##provider: GENCODE
    ##contact: gencode-help@sanger.ac.uk
    ##format: gtf
    ##date: 2016-07-15
    chr1    HAVANA  gene    11869   14409   .   +   .   gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";
    chr1    HAVANA  transcript  11869   14409   .   +   .   gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_status "KNOWN"; transcript_name "DDX11L1-002"; level 2; transcript_support_level "1"; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";

## Data repositories

Raw sequencing data (FASTQ files) and sample metadata:

- ENA (https://www.ebi.ac.uk/ena/browser/guides). 
- NCBI (https://www.ncbi.nlm.nih.gov/) (GEO/SRA).

Genome assemblies and annotation:

- ENSEMBL (https://www.ensembl.org/index.html)
- NCBI

## High Performance Computing (HPC)

For login to the CCC open a Linux terminal use SSH conexion.

    - accesing only via terminal
        
        ssh user@login1.ccc.uam.es

    - accesing with graphical interface support

        ssh -X user@login3.ccc.uam.es

Software on the CCC is already installed and managed through the module system.

You can use the "module" command.

```bash
module avail                      # lists all software available
module load software_name         # load a specific software into your environment
module list                       # displays all software currently loaded
module unload software_name       # unloads a specific sofware module
module purge                      # unloads all currently loaded software
```
At the CCC you have access to the following directories:

```bash
/home/user                        # your home directory 
/scratch/user                     # temporary high-performance storage
/home/projects/your_project       # shared project directory
```

## Visual Studio Code - SSH extension

Powerful code editor that supports many programming languages such as Bash, Python, R, and others. 

To work directly on the CCC from your local machine, you can use the remote SSH extension.

**VSCode configuration for command line only (no graphics)**

*Configuration parameters*:

    Host "alias_host"
        HostName login1.ccc.uam.es
        User user_name

**VSCode configuration with graphic support (X11 Forwarding)** \
Note: In Windows, this configuration is not enough to see graphics.

*Configuration parameters*:

    Host "alias_host"
        HostName login3.ccc.uam.es
        User user_name
        ForwardX11 yes

## PuTTY with X11 Forwarding (PuTTY + XLaunch)

For Windows, you can connect to a remote server using PuTTY with X11 forwarding, using VcXsrv (XLaunch).

Download VcXsrv from:

    https://sourceforge.net/projects/vcxsrv/

Install VcXsrv and launch XLaunch with the default settings. This will start the X server in the background (you should see its icon in the system tray).

Download PuTTY from:

    https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html

Install and open PuTTY. Configure a session with the following settings: 

    Host: login3.ccc.uam.es
    Port: 22
    Connection type: SSH

In the left menu, go to:

    Connection > SSH > X11

Check "Enable X11 forwarding"

Make sure XLaunch is runing before opening PuTTY, otherwise X11 applications will not display.

You can save the session setting by entering a name in the "Saved Sessions" and click "Save". You can then load this session for future connections.

## Conda environments

Conda is an environment manager that allow to create isolated environments and avoid conflicts between packages dependecies
Also to install specific versions of software and libraries.
It is included in distribution such as Anaconda, Miniconda and Bioconda.

On the CCC, Conda is available through the Miniconda module.

    module load miniconda
    
Before using Conda, you must initialize it:

    conda init
    source ~/.bashrc

Now you can work with Conda. These are some usefull commands:

```bash
conda env list                          # shows available environments
conda create --name new_env             # creates an environment
conda env create -f environment.yml     # creates an environment using a config file
conda activate env_name                 # activates an environment
conda install software_name             # install package in the active environment
conda list                              # displays all packages installed
conda deactivate                        # returns to base environment
conda export > environment.yml          # export your current environment 
```



# RNASEQ ANALYSIS 

## 1. Working environment setup

```
# Log in the CCC allowing graphics

    ssh -X user@login3.ccc.uam.es

# Move to the scratch dir

    cd /scratch/user

# Create a new dir for your analysis project

    mkdir my_rna_project

# Initiate conda 

    module load miniconda
    conda init
    source ~/.bashrc

# Check available conda environments 

    conda env list

# Import conda environments
# for rna-seq analysis

    wget "https://path/to/rnaseq_env.yml"

# for differential expression and GSEA analysis

    wget "https://path/to/deseq2_env.yml"

# Create conda environments from the .yml files

    conda env create -f rnaseq_env.yml
    conda env create -f deseq2_env.yml

# Activate an environment

    conda activate rnaseq_env

```

## 2. Sequencing data download

Go to ENA browser (https://www.ebi.ac.uk/ena/browser/home).

Enter the study name (PRJNA998264) and get the path to the reads for the samples of interest (SRR25411155, SRR25411161, SRR25411166, SRR25411167). 

Select those samples and click "Get download script"

```
# Make sure you are on the project folder

    pwd         # should be /scratch/user/my_rna_project

# Create a new folder

    mkdir -p 01_reads/01_raw_reads
    cd 01_reads/01_raw_reads

# Download sequencing read files (FASTQ) for your samples
    
    wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/055/SRR25411155/SRR25411155.fastq.gz # -nc, avoids overwritting files taht already exits
    wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/067/SRR25411167/SRR25411167.fastq.gz
    wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/066/SRR25411166/SRR25411166.fastq.gz
    wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/061/SRR25411161/SRR25411161.fastq.gz

# Decompress the files

    gzip -d -v *.fastq.gz

# Examine one of the files to check the format

    head SRR25411155.fastq

# Count the number of lines to get the number of reads

    wc -l SRR25411155.fastq

# Subsample the reads for the training
# (Optional) Check how the command will look like before running it using "echo"
#   
#   for file in *.fastq
#   do
#       echo "seqtk sample -s100 ${file} 0.1 > ${file%%.fastq}_subsampled.fastq"
#   done

    for file in *.fastq
    do
        seqtk sample -s100 ${file} 0.1 > ${file%%.fastq}_subsampled.fastq
    done

```

## 3. Reference genome preparation
Go to ENSEMBL (https://www.ensembl.org/index.html)

Look for the organism of interest and get the link from "Download FASTA" and "Download GTF".

```
# Create a folder to store the reference genome data

    mkdir /scratch/user/my_rna_project/02_reference
    cd /scratch/user/my_rna_project/02_reference

# Download the reference genome sequence (FASTA) and annotation (GTF) files.

    wget https://ftp.ensembl.org/pub/release-115/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.gz
    wget https://ftp.ensembl.org/pub/release-115/gtf/danio_rerio/Danio_rerio.GRCz11.115.gtf.gz

# Decompress the files

    gzip -d -v *.gz

# Build a HISAT2 index for the genome

    module load hisat2/2.2.1
    hisat2-build -p 24 Danio_rerio.GRCz11.dna_sm.primary_assembly.fa Drerio

```

## 4. Reads QC assessment
Assess the sequencing quality of the downloaded sequences.

```
# Create a dir to store the QC reports

    mkdir /scratch/my_rna_project/01_reads/01_raw_fastqc
    cd /scratch/my_rna_project/01_reads/01_raw_fastqc

# Run FASTQC

    module load fastqc
    fastqc *_subsampled.fastq -t 6 --outdir /home/proyectos/cbmngs/cursoRNA202601/my_rna_project/01_reads/01_raw_fastqc/

# Run multiQC to aggregate FASTQC reports

    multiqc .

# Inspect quality reports

    firefox *.html

```

## 5. Reads trimming
Trimme the reads to remove adapter sequences and to get rid of very short and low quality reads.

```
# Create folers to store the trimmed reads and the QC reports after trimming

    mkdir /scratch/user/my_rna_project/01_reads/02_trimmed_reads
    mkdir /scratch/user/my_rna_project/01_reads/02_trimmed_fastqc

    cd /scratch/user/my_rna_project/01_reads/01_raw_reads

# Run trim_galore

    for file in *_subsampled.fastq
    do 
	    trim_galore --length 20 ${file} --fastqc --fastqc_args "--outdir /scratch/user/my_rna_project/01_reads/02_trimmed_fastqc" -j 6 -o /scratch/user/my_rna_project/01_reads/02_trimmed_reads/
    done

# Check the FASTQC reports for trimmed reads
    
    cd ../02_trimmed_fastqc/
    multiqc .
    firefox multiqc_report.html

# Compare the fastqc report for the same sample before and after trimming

    firefox SRR25411155_subsampled_trimmed_fastqc.html ../01_raw_fastqc/SRR25411155_subsampled_fastqc.html 

```

## 6. Splice aware alignment (HISAT2)

```
# Create a dir to store alignment files

    mkdir /scratch/user/my_rna_project/03_hisat2

# Run HISAT2 on trimmed reads. Note: if you don't know if your RNAseq experiment is stranded or not, there are softwares for checking it.

    module load hisat2/2.2.1
    module load samtools/1.9

    cd /scratch/user/my_rna_project/01_reads/02_trimmed_reads

    for file in *.fq
    do
	    sample_name=${file%%_subsampled_trimmed.fq}
	    echo "Processing sample ${sample_name}"
	    hisat2 -x /scratch/user/my_rna_project/02_reference/Drerio -q -k 1 -p 24 -U ${file} -S /scratch/user/my_rna_project/03_hisat2/${sample_name}.sam 2>> /scratch/user/my_rna_project/03_hisat2/${sample_name}_hisat2_summarymetrics.txt
    done

# Move to the alignment folder

    cd /scratch/user/my_rna_project/03_hisat2

# Transform the SAM file in compressed BAM file format and sort them by position
    
    for file in *.sam
    do
	    sample_name=${file%%.sam}
	    samtools view -S -b ${sample_name}.sam --threads 24 > ${sample_name}.bam
	    samtools sort ${sample_name}.bam -o ${sample_name}_sorted.bam -@ 24
	    samtools index ${sample_name}_sorted.bam -@ 24
    done

# Check your alignments in IGV

    module load igv
    igv.sh

```

## 7. Reads per gene quatification (featureCounts)

```

# Create a dir to store quantification results

    mkdir /scratch/user/my_rna_project/04_featureCounts

# Run featureCounts

    module load featureCounts/2.1.1 
    cd /scratch/user/my_rna_project/03_hisat2

    for file in *_sorted.bam
    do
	    sample_name=${file%%_sorted.bam}
	    featureCounts -T 5 -t exon -g gene_id -a /scratch/user/my_rna_project/02_reference/Danio_rerio.GRCz11.115.gtf -o /scratch/user/my_rna_project/04_featureCounts/${sample_name}_ftC.tsv ${file} 
    done

```

## 8. Differential expression (DESEQ2)

```
# Set up a new conda environment

    conda deactivate
    module purge # to remove all previously loaded modules
    conda activate deseq2_env

# Create a folder to run the analysis and store the results

    mkdir /scratch/user/my_rna_project/05_deseq2

# Format featurecounts output file so they can be used as input in DESEQ2

    cd /scratch/user/my_rna_project/04_featureCounts

    for file in *_ftC.tsv
    do
	    sample_name=${file%%_ftC.tsv}
	    awk -F'\t' '{print $1 "\t" $7}' "${file}" | tail -n+3 > "/scratch/user/my_rna_project/05_deseq2/${sample_name}.tsv"
    done

# Make a configuration TSV file for DESEQ2.
# it should have three columns: sample  file    condition
# samples with the baseline condition (e.g. WT) should be placed on the first rows.

# Run deseq2

    cd /scratch/user/my_rna_project/05_deseq2
    Rscript --vanilla deseq2_script.R -c deseq2_config.tsv -d org.Dr.eg.db

# Visualize output .tiff files

    xdg-open deseq2_pca_labels.tiff

# Visualize interactive .html files

    firefox deseq2_volcano_mutant_vs_WT.html

```

## 9. Pathway enrichment analysis (GSEA with GO terms)

```
# Create a dir to store pathway enrichment results

    mkdir /scratch/user/my_rna_project/06_GSEA

# Filter DESEQ2 differential expression output

    tail -n+2 deseq2_all_mutant_vs_WT.tsv | cut -f1,7 | grep -P -v "\tNA" > /scratch/user/my_rna_project/06_GSEA/mutant_vs_WT_for_gsea.tsv

# Run GSEA 

    Rscript --vanilla gseaGO_analysis.R -i mutant_vs_WT_for_gsea.tsv -p Mutant_WT -o org.Dr.eg.db

```
