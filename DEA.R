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
