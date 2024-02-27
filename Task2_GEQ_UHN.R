# Title: Gene expression quantification on RNA Seq Fastq files
# File: Task2_GEQ_UHN
# Project: UHN_Assignment

# Working Dir: setwd("C:/Users/surab/Desktop/UHN_Tasks")

# Install necessary packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("ShortRead", "Rsamtools", "GenomicFeatures", "DESeq2"))

# Load required libraries
library(ShortRead)
library(Rsamtools)
library(GenomicFeatures)
library(DESeq2)

# Read FASTQ files
fastq_files <- c("1A_1.fastq", "1A_2.fastq", "1B_1.fastq", "1B_2.fastq", "2A_1.fastq", "2A_2.fastq", "2B_1.fastq", "2B_2.fastq", "2C_1.fastq", "2C_2.fastq")
reads <- readFastq(fastq_files)

# Perform quality control (optional)
# You can run FastQC from within RStudio using system commands

# Align reads to a reference genome
# Example with Bowtie2
# system("bowtie2-build reference_genome.fa reference_genome")
# system("bowtie2 -x reference_genome -U sample1.fastq -S sample1.sam")
# (Replace 'reference_genome' with your actual reference genome)

# Convert SAM files to BAM files
# samtools <- system.file("exec", "samtools", package="Rsamtools")
# system(paste(samtools, "view -bS", "sample1.sam", "-o", "sample1.bam"))

# Load gene annotation (GTF file)
gtf_file <- "annotation.gtf"
genes <- makeTxDbFromGFF(gtf_file)

# Count reads mapping to genes
# Replace "sample1.bam" with actual aligned BAM file
bamfile <- BamFile("sample1.bam", yieldSize=100000)
txdb <- genes
tx_by_gene <- exonsBy(txdb, by="gene")
counts <- summarizeOverlaps(tx_by_gene, bamfile, mode="Union", singleEnd=TRUE,
                            ignore.strand=TRUE, fragments=TRUE)
counts <- assays(counts)$counts

# Create a DESeqDataSet object
colData <- data.frame(condition = c("condition1", "condition2"))
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ condition)

# Perform differential expression analysis
dds <- DESeq(dds)
res <- results(dds)


save.image(file = "RNASeq_GEQ_UHN.RData")
load("RNASeq_GEQ_UHN.RData")
