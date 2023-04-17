###########################################################
#                                                         #
# Differential gene expression analysis (RNA-seq) with R  #
#                                                         #
#                                                         #
#                                                         #
###########################################################

# Análisis de genomas y transcriptomas con plataforma NGS
# Master Biotecnología
# APR 2023

# Jose V. Die, Dep Genetics, UCO
# jose.die@uco.es


# Dependencies
library("Rbowtie2")
library("Rsamtools")
library("Rsubread")
library("DESeq2")


# Organize the Project -------------------------------------------------------------------------------
# Download genome data
download.file("https://www.dropbox.com/s/4ft480eky7kghzw/Saccharomyces_cerevisiae_genome.fa.gz?dl=1", 
              destfile = "Saccharomyces_cerevisiae_genome.fa.gz")

download.file("https://www.dropbox.com/s/qtaret1hrbvw2xb/Saccharomyces_cerevisiae_genome.gff3.gz?dl=1", 
              destfile = "Saccharomyces_cerevisiae_genome.gff3.gz")

# Download sequencing reads
download.file("https://www.dropbox.com/s/v06um7vt9ojdf42/NGS_MB_SRR9336468_1.fastq.gz?dl=1", 
              destfile = "1M_SRR9336468_1.fastq.gz")

download.file("https://www.dropbox.com/s/kgxfth8ra675ccu/NGS_MB_SRR9336468_2.fastq.gz?dl=1", 
              destfile = "1M_SRR9336468_2.fastq.gz")

# Uncompress genome 
dir.create("genes")
system("time gunzip -k Saccharomyces_cerevisiae_genome.fa.gz")
system("time gunzip -k Saccharomyces_cerevisiae_genome.gff3.gz")
system("mv Saccharomyces_cerevisiae_genome.* genes")

# Uncompress reads
dir.create("fastq")
system("time gunzip -k *.fastq.gz")
system("mv *.fastq *.gz fastq/")


# STEP 1: Create genome index for bowtie2 with "bowtie2_build" ---------------------------------------------
bowtie2_build(references = "genes/Saccharomyces_cerevisiae_genome.fa",
              bt2Index = "genesScerevisiae_genome",
              overwrite = TRUE)

# Index the genome fasta file with Rsamtools
indexFa("genes/Saccharomyces_cerevisiae_genome.fa") # crea el .fai



# STEP 2: Align FASTQ files against indexed genome with "bowtie2" --------------------------------------------
dir.create("bam")

bowtie2(bt2Index = "genes/Sc_genome", 
        samOutput = "bam/1M68_pH5_0.04CO2_R1.sam", 
        seq1 = "fastq/1M_SRR9336468_1.fastq", 
        seq2 = "fastq/1M_SRR9336468_2.fastq", "--threads=3")       



# STEP 3 : Convert SAM files into BAM (and indexes .bai) with "asBam" --------------------------------------------
asBam(file = "bam/1M68_pH5_0.04CO2_R1.sam", destination = "bam/basal_1")

# make some tidy : una vez se tiene el BAM, el fichero(s) SAM se puede borrar porque son los mas prescindibles y los que mas ocupan
file.remove("SRR9336468.sam")  

# download the rest of BAM files : 
#https://www.dropbox.com/s/pr3vqfj9r2c35xg/bam_files.zip?dl=1



# STEP 4: Create metadata file for the Saccharomyces count matrix --------------------------------------------
# Create sample vector
sampleID = paste0("SRR93364", 68:76)

# Create condition vector
conditions = rep(c("basal", "indust1", "indust2"), each = 3)

# Create description vector
description = rep(c("pH5_0.04CO2", "pH5_5CO2", "pH3_0.04CO2"), each = 3)

# Create data frame
saccha_metadata <- data.frame(sampleID, conditions, description)

# Assign the row names of the data frame
rownames(saccha_metadata) = paste(conditions, 1:3, sep = "_")

# Save the df into a csv file 
dir.create("data")
write.csv(saccha_metadata, file = "data/saccha_metadata.csv")



# STEP 5: Quantify "genes" with "featureCounts" -------------------------------------------------------------
# store all paths of bam files into "bamFiles" R object (use "list.files"
list.files(path = "bam", pattern = "*.bam$", full.names = TRUE)
bamFiles <- list.files(path="bam", pattern="*.bam$", full.names=TRUE)

# quantify 'genes'
data <- featureCounts(
  files           = bamFiles,
  annot.ext       = "genes/Saccharomyces_cerevisiae_genome.gff3",
  isGTFAnnotation = TRUE,
  GTF.featureType = "gene",
  GTF.attrType    = "ID",
  isPairedEnd     = TRUE,
  requireBothEndsMapped = TRUE,
  minMQS          = 30,
  strandSpecific  = 2,
  nthreads = 4) 


# store the results in the object `raw_counts``
dir.create("results")
write.table(data$counts, file = "results/raw_counts" , sep = "\t")

# .tsv means "tab separated values"
#raw_counts = read.table("results/raw_counts")



#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#------- CONTINUAR desde AQUI : 

# trabajar con la clase los siguientes conceptos : 
# mapping gen YAL063C in the last 4 samples 
raw_counts["gene:YAL063C", 6:9]

# does the mapping match what I see in IGV ?


# open IGV : https://igv.org/
# load Saccharomices genome
# search : YAL063C (v. desktop no permite busqueda ; v. app si lo encuentra y 
# reconoce que el nombre biologico del gen YAL063C es FLO9)
# v. desktop : buscar por : chrI:16,439-36,407

# load samples (bam) 73-74-75-76 (f. ex.)
# color by first-pair-strand
# set track height : 200



# # STEP 9: Prepare data for the DESeq2 workflow
# -----------------------------------------------------------------------------  

# Organize the data   
# DESeq2 requires the sample names in the metadata and counts dataset to be in the 
# same order. Therefore, the row names in the metadata need to be in the same order 
# as the column names of the count data (= mismo orden que los ficheros bam).

# - raw counts of the number of reads aligning each gene 
# - associated sample metadata 

raw_counts = read.table("results/raw_counts")
samples  <- read.csv("data/saccha_metadata.csv", row.names = 1)

colnames(raw_counts)
rownames(samples)

all(paste0(rownames(samples), ".bam") == colnames(raw_counts))
identical(paste0(rownames(samples), ".bam") , colnames(raw_counts))


# STEP 10: Use "DESeqDataSetFromMatrix" to create a DESeq2 object
# -----------------------------------------------------------------------------  
dds <- DESeqDataSetFromMatrix(countData = raw_counts, 
                              colData = samples, 
                              design = ~ conditions)

save(dds, file = "results/dds.rda")
#load(file = "results/dds.rda")

#methods(class = "DESeqDataSet")
#levels(dds$conditions)



# STEP 11: Normalize raw counts for library size
# -----------------------------------------------------------------------------  

# The raw counts for each sample are divided by the associated sample-specific 
# size factor for normalization. To view the size factors used for normalization, 
# we can use the sizeFactors() function.
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
# colSums(counts(dds))
# 
# plot(sizeFactors(dds),colSums(counts(dds)))
# abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))
# 
# rs <- rowSums(counts(dds))
# boxplot(log2(counts(dds)[rs > 0,] + 1)) # not normalized
# log.norm.counts <- log2(counts(dds, normalized=TRUE) + 1)
# boxplot(log.norm.counts[rs > 0,]) # normalized

#scatterplot of log normalized counts against each other.
# Note the fanning out of the points in the lower left corner, for points less 
#than 2^5=32.
plot(log.norm.counts[,1:2], cex=.1)

# Once the size factors have been calculated and added to the DESeq2 object, 
# the normalized counts can be extracted. 
# If the default was left as normalized = FALSE, then we would extract the raw 
# counts from the object.
norm_counts <- counts(dds, normalized = TRUE)
#boxplot(log2(norm_counts+1))
write.csv(norm_counts, file = "results/norm_counts.csv")



# STEP 12: Unsupervised Clustering Analyses
# -----------------------------------------------------------------------------  

# Transform the normalized counts : VST-transformed normalized counts
vst_counts <- vst(dds, blind = TRUE)

# Extract the matrix of transformed counts
vsd <- assay(vst_counts)

# Compute the correlation values between samples
vsd_cor <- cor(vsd) 
#View(vsd_cor)


# Hierarchical analysis 
#BiocManager::install("pheatmap")
#samples  <- read.csv("data/saccha_metadata.csv", row.names = 1)
pheatmap::pheatmap(vsd_cor, annotation = dplyr::select(samples, conditions))

# hemos hecho subset de 1M
# si tuvieramos las lecturas completas, puede que se hubieran agrupado bien 
# hay que decidir si se deja fuera una muestra, si se repite, ...

# ex. : simulate other metadata such as the sample collection date
#samples$collection <- c(rep("April", 2), rep("June", 2), rep("April", 5))
#pheatmap::pheatmap(vsd_cor, annotation = dplyr::select(samples, collection))


# Plot PCA
plotPCA(vst_counts, intgroup = "conditions")

#This is great since it seems that a lot of variation in gene expression in the 
#data set can likely be explained by the differences between sample groups. 



# STEP 13: Run the DESeq2 analysis : model fitting for gene expression data
# -----------------------------------------------------------------------------  
dds <- DESeq(dds)

#Now that we have run the DE analysis, we could explore our results. 
#However, before proceeding, we should explore how well our data fit the model.
#datacamp : deseq2model video

# Plot dispersions
plotDispEsts(dds)



# STEP 14: Run the DESeq2 model contrast to extract the results
# -----------------------------------------------------------------------------  
P5C5vsP5C004 <- results(dds, 
                         contrast=c("conditions", "indust1", "basal"), 
                         alpha = 0.05)


plotMA(P5C5vsP5C004, colSig = "red")

#save(P5C5vsP5C004, file = "results/P5C5vsP5C004.rda")
#load("results/P5C5vsP5C004.rda")


write.csv(as.data.frame(P5C5vsP5C004), 
          file="results/P5C5vsP5C004_DEG.csv")

write.table(P5C5vsP5C004, file="results/P5C5vsP5C004_DEG.tsv", sep="\t", row.names=FALSE)



# STEP 15: DESeq2 results table 
# -----------------------------------------------------------------------------  
mcols(P5C5vsP5C004)
summary(P5C5vsP5C004)
head(P5C5vsP5C004)


# set a log2 threshold to extract genes with biological relevance 
P5C5vsP5C004 <- results(dds, 
                        contrast=c("conditions", "indust1", "basal"), 
                        alpha = 0.05, 
                        lfcThreshold = 0.32)     # 1.25-fold-change


# Extract the significant genes (p-adjusted < 0.05) : 'subset'
dim(P5C5vsP5C004)
dim(subset(P5C5vsP5C004, padj < 0.05))

dim(subset(P5C5vsP5C004, log2FoldChange > 0 & padj < 0.05))
dim(subset(P5C5vsP5C004, log2FoldChange < 0 & padj < 0.05))

subset(P5C5vsP5C004, log2FoldChange > 0 & padj < 0.05)
subset(P5C5vsP5C004, log2FoldChange < 0 & padj < 0.05)




# STEP 16: Visualizing results
# -----------------------------------------------------------------------------  
# Subset normalized counts to significant genes
norm_counts <- read.csv(file = "results/norm_counts.csv", row.names = 1)

degs <- rownames(subset(P5C5vsP5C004, padj < 0.05))
sig_norm_counts <- norm_counts[degs, ] # tabla final para heatmap


# Choose a color palette from RColorBrewer
library(RColorBrewer)
heat_colors <- brewer.pal(6, "YlOrRd")


# Visualization 1 : MA plot
plotMA(P5C5vsP5C004, colSig = "red")


# Visualization 2 : heatmap
#samples  <- read.csv("data/saccha_metadata.csv", row.names = 1)
pheatmap::pheatmap(sig_norm_counts, 
                   color = heat_colors, 
                   cluster_rows = TRUE, 
                   show_rownames = FALSE, 
                   annotation = dplyr::select(samples, conditions), 
                   scale = "row")


# Visualization 2 : volcano plot

CONTINUAR


# STEP 17 : Annotation of DEGs
# -----------------------------------------------------------------------------  

# Download genome annotation 
download.file("https://www.dropbox.com/s/clp3eson00mlca3/sacCer_annotation.csv?dl=1", 
              destfile = "data/Saccharomyces_annotation.csv")

# Read the annotation
annot <- read.csv("data/Saccharomyces_annotation.csv")


# Table with DEGs
res <- as.data.frame(subset(P5C5vsP5C004, padj < 0.05))

# Load function
library(dplyr)
final_table <- function(df) {
  
  genes = sub("gene:", "", rownames(df))
  
  tmp <- df %>% 
    mutate('Locus.tag' = genes) %>% 
    select('Locus.tag', everything()) %>% 
    as_tibble() %>% 
    inner_join(annot, by = 'Locus.tag') %>% 
    rename(gene = Locus.tag, chr = X.Name) %>% 
    select(-GeneID) %>% 
    select(gene, Protein.Name, everything())
  
  tmp
}


# Call function `final_table` and get final table DEGS with annotation 
res_final <- final_table(res)

res_final


res_final %>% arrange(desc(log2FoldChange ))











# Upload text-tabulated file into FIESTA viewer to create an interactive MA plot:
#
# https://bioinfogp.cnb.csic.es/tools/FIESTA/index_server.php
#
###############################################
#                                             #
# Gene Identifier     : GENEID                #
#                                             #
# X axis (horizontal) : log2BaseMean          #
#                                             #
# Y axis (vertical)   : log2Ratio             #
#                                             #
# Error Y (vertical)  : STDERR_log2Ratio      #
#                                             #
# Spot Colors         : MA plot               #
#                                             #
# Press [Next]                                #
#                                             #
###############################################

# Wait few seconds and download the ".zip" file with all results.
# Uncompress ".zip" folder in your hard disk an open "index.html" with any browser.
#
# Filter (eg. "log2Ratio>0" and "padjust < 0.05"),
# sort results by padjust and save the table
#
# PROPOSED TASKS: 
# ==============
# 1. Filter and save down-regulated genes in comparison P5C50vsP5C004
#
# 2. Get up-regulated an down-regulated genes for comparison: "pH=3 CO2=0.04%" versus "pH=5 CO2=0.04%"
# P3C004vsP5C004 <- results(dds, contrast=c("conditions", "p3c0.04", "p5c0.04"))

# 3. Play around with different filtering thresholds (padjust < 0.1; padjust <0.01; others...)
#
library(tibble)
library(dplyr)
as_tibble(P3C004vsP5C004) %>% 
  filter(padj < 0.1)

as_tibble(P3C004vsP5C004) %>% 
  filter(padj < 0.01)




# 4. (DIFFICULT): Create a text-tabulated file with raw counts. Tip: assays(dds)$counts
#
# 5. (DIFFICULT): Create a text-tabulated file with normalized counts. Tip: assays(dds)$mu
#           NOTE: It is illustrative to look at normalized counts of differentially
#               expressed genes with good pvalues.
#
# 6. (DIFFICULT): add columns with normalized counts to the FINAL table (create file).
#
# Useful links:
# =============
# Illumina Sequencing Technology (Many thanks to Victoria Castro Illana!)
# https://www.youtube.com/watch?v=womKfikWlxM&frags=pl%2Cwn
#
# All about DESeq2 bioconductor package:
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#
# About MA plots (wikipedia):
# https://en.wikipedia.org/wiki/MA_plot
#
# padjust (FDR calculation) clearly explained!!
# https://statquest.org/tag/fdr/


########################################################################
#
# strandSpecific = 0 means "no strand specific"
# strandSpecific = 1 means "fragment _1 is direct strand"
# strandSpecific = 2 means "fragment _1 is reverse strand"
#
# mRNA: 5'-AUGCGUCGACGUUUAGCAUGUCGAUGCUGAUUGCGAUCGAUGCCGACAAGAUAGUAA-3'
# 
#              (retrotranscription and second strand generation)
# 
# cDNA: 5'-ATGCGTCGACGTTTAGCATGTCGATGCTGATTGCGATCGATGCCGACAAGATAGTAA-3'
#          |||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#       3'-tacgcagctgcaaatcgtacagctacgactaacgctagctacggctgttctatcatt-5'
#                                                    <-------------- 
#                                                  (note the polarity)     
#
#
# strandSpecific=0:  reads_?.fastq:              reads_?.fastq:
#                    5'-ATGCGTCGACGTTTA-3'       5'-ttactatcttgtcgg-3'                   
#                                                   -------------->
#                                                                                                                                                     
# strandSpecific=1:  reads_1.fastq:              reads_2.fastq:
#                    5'-ATGCGTCGACGTTTA-3'       5'-ttactatcttgtcgg-3'                   
#                                                   -------------->
#                                                                                                                     
# strandSpecific=2:  reads_1.fastq:              reads_2.fastq:
#                    5'-ttactatcttgtcgg-3'       5'-ATGCGTCGACGTTTA-3'                   
#                       -------------->
#
########################################################################


