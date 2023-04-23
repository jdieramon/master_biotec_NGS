###########################################################
#                                                         #
#    Mapping  Reads                                       #
#                                                         #
#                                                         #
###########################################################


# Mappig Reads to a Genome 
# Análisis de genomas y transcriptomas con plataforma NGS
# Master Biotecnologia
# APR 2023

# Jose V. Die, Dept. Genetics, UCO
# jose.die@uco.es


# Paired-end reads
#	FFFF_1.fastq.gz
#	FFFF_2.fastq.gz

# .FASTQ (.FQ):       short-reads (eg. illumina)
# .FASTA (.FA,.FNA):  biological sequences (eg. chromosomes)
# .GFF (.GFF3, .GTF): genomic elements coordinates (eg. genes)

# Experiment:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133136

# Dependencies
library(Rbowtie2)
library(Rsamtools)


# Organize the Project
# dir.create("data")    #(To create a directory/folder in the indicated path)
# dir("data")           #(To see files and folders in the indicated path)



# Download data ----------------------------------------------------------------
file.exists("1M_SRR9336468_1.fastq.gz")
file.exists("1M_SRR9336468_2.fastq.gz")
file.exists("Saccharomyces_cerevisiae_genome.gff3.gz")
file.exists("Saccharomyces_cerevisiae_genome.fa.gz")


download.file("https://www.dropbox.com/s/v06um7vt9ojdf42/NGS_MB_SRR9336468_1.fastq.gz?dl=1", 
              destfile = "data/1M_SRR9336468_1.fastq.gz")

download.file("https://www.dropbox.com/s/kgxfth8ra675ccu/NGS_MB_SRR9336468_2.fastq.gz?dl=1", 
              destfile = "data/1M_SRR9336468_2.fastq.gz")

download.file("https://www.dropbox.com/s/qtaret1hrbvw2xb/Saccharomyces_cerevisiae_genome.gff3.gz?dl=1", 
              destfile = "data/Saccharomyces_cerevisiae_genome.gff3.gz")

download.file("https://www.dropbox.com/s/4ft480eky7kghzw/Saccharomyces_cerevisiae_genome.fa.gz?dl=1", 
              destfile = "data/Saccharomyces_cerevisiae_genome.fa.gz")



# Check files
system("ls -lh")
dir()
list.files()


# Uncompress genome fastq gz files with "gunzip"
#R.utils::gunzip("Saccharomyces_cerevisiae_genome.gff3.gz", remove = FALSE)
system("time gunzip -k data/Saccharomyces_cerevisiae_genome.fa.gz")
system("time gunzip -k data/Saccharomyces_cerevisiae_genome.gff3.gz")

# Uncompress reads fastq gz files with "gunzip"
system("time gunzip -k data/NGS*.fastq.gz")

# Uncompress all fastq gz files with "gunzip"
system("time gunzip -k data/*.fastq.gz")


# Organize the Project
# dir.create("data/chromosomes")
# dir.create("data/fastq")
system("mkdir data/chromosomes")
system("mkdir data/fastq")
system("mkdir results")
system("mv data/NGS*.fastq data/fastq/")
system("mv data/Saccharomyces_cerevisiae_genome.fa data/chromosomes/")
system("mv data/Saccharomyces_cerevisiae_genome.gff3 data/chromosomes/")
#system("rm *.gz")

system("ls -lh")
system("ls -lh data/chromosomes")
system("ls -lh data/fastq")




## GFF3 & FASTQ files ----------------------------------------------------------
# Look at the genomic annotation file (gff3)
file.show("data/chromosomes/Saccharomyces_cerevisiae_genome.gff3")
system("less -S data/chromosomes/Saccharomyces_cerevisiae_genome.gff3")

# Look at one of the two *.fastq files
system("wc -l data/fastq/1M_SRR9336468_1.fastq")
system("less -S fastq/1M_SRR9336468_1.fastq")

# Look at the two reads from the same cDNA fragment 
system("head -n 4 data/fastq/1M_SRR9336468_1.fastq")
system("head -n 4 data/fastq/1M_SRR9336468_2.fastq")



## BLAST genome ------------------------------------------------------------
# blast genome : saccharomices
# https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&USER_FORMAT_DEFAULTS=on&SET_SAVED_SEARCH=true&PAGE=MegaBlast&PROGRAM=blastn&GAPCOSTS=0%200&MATCH_SCORES=1,-2&BLAST_SPEC=Assembly&DATABASE=genomic/559292/GCF_000146045.2&BLAST_PROGRAMS=megaBlast&MAX_NUM_SEQ=100&SHORT_QUERY_ADJUST=on&EXPECT=0.05&WORD_SIZE=28&REPEATS=4932&TEMPLATE_TYPE=0&TEMPLATE_LENGTH=0&FILTER=L&FILTER=R&FILTER=m&EQ_MENU=Enter%20organism%20name%20or%20id--completions%20will%20be%20suggested&PROG_DEFAULTS=on&SHOW_OVERVIEW=on&SHOW_LINKOUT=on&ALIGNMENT_VIEW=Pairwise&MASK_CHAR=2&MASK_COLOR=1&GET_SEQUENCE=on&NUM_OVERVIEW=100&DESCRIPTIONS=100&ALIGNMENTS=100&FORMAT_OBJECT=Alignment&FORMAT_TYPE=HTML



## Genome alignment ------------------------------------------------------------
#bowtie2 paper : https://www.nature.com/articles/nmeth.1923
#bowtie2 software homepage : http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
#BiocManager::install("Rbowtie2")
#library(Rbowtie2)

# STEP 1: Create genome index for bowtie2 with "bowtie2_build":
bowtie2_build(references = "data/chromosomes/Saccharomyces_cerevisiae_genome.fa",
              bt2Index = "data/chromosomes/Scerevisiae_genome",
              overwrite = TRUE)

# Index the genome fasta file with Rsamtools
indexFa("data/chromosomes/Saccharomyces_cerevisiae_genome.fa") # crea el .fai


# STEP 2: Align FASTQ files against indexed genome with "bowtie2":
bowtie2(bt2Index = "data/chromosomes/Scerevisiae_genome",
        samOutput = "SRR9336468.sam",
        seq1 = "data/fastq/1M_SRR9336468_1.fastq",
        seq2 = "data/fastq/1M_SRR9336468_2.fastq",
        overwrite = TRUE,
        "--threads=3")


# STEP 3 : Convert SAM files into BAM (and indexes .bai) with "asBam"
#BiocManager::install("Rsamtools")
#library(Rsamtools)
asBam("SRR9336468.sam")


# make some tidy
# una vez se tiene el BAM se puede borrar el SAM porque son los mas prescindibles y los que mas ocupan
file.remove("SRR9336468.sam")                  





## IGV  ------------------------------------------------------------
# IGV software homepage : https://www.broadinstitute.org/igv/
# https://igv.org/
  
# Load genome : Genome / Local File / chromosomes /  .fa + .fai
# Load genes : Tracks / Local File (gff3)
#View options : 
  #Set track height / 100
  #Expanded / Collapse
# Load reads : Tracks / Load File / select bam + index (bai)



## Filtering reads by mapping quality
# ------------------------------------------------------------------------------
bamfile <- BamFile("SRR9336468.bam")
bamfile
seqinfo(bamfile)

# set the filter condition
# build the filtered : .bam + .bai files 
param = ScanBamParam(mapqFilter = 35) # min=35 * 
dest = filterBam(file = bamfile, destination = "results/mapq.bam", param = param) 

#dest <- BamFile("results/mapq.bam")
countBam(bamfile)
countBam(dest)

#scanBam generates a list (length = 1) with 13 elements
aln <- scanBam(bamfile)
aln35 = scanBam(dest)


#names(aln[[1]])
aln[[1]]$seq[1:10]
aln[[1]]$mapq[1:10]
summary(aln[[1]]$mapq)

summary(aln[[1]]$mapq)
summary(aln35[[1]]$mapq)

sum(aln[[1]]$mapq < 35, na.rm = TRUE)
sum(aln35[[1]]$mapq < 35, na.rm = TRUE)

