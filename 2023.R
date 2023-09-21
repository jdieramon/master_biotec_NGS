
###########################################################
#                                                         #
#  Sequence reads assessment
#                                                         #
#                                                         #
#                                                         #
###########################################################

# Analisis de genomas y transcriptomas con plataforma NGS
# Master Biotecnologia
# APR 2022

# Jose V. Die, Dep Genetics, UCO
# jose.die@uco.es


#BiocManager::install("ShortRead")

# Dependencies
library(ShortRead)



# Download files 
download.file("https://www.dropbox.com/s/6un05pdvx2join5/Saccharomyces_cerevisiae_genome.fa.gz?dl=1", 
              destfile = "Saccharomyces_cerevisiae_genome.fa.gz")


download.file(url = "https://www.dropbox.com/s/svch71j19zfhntm/SRR11397715.fastq.gz?dl=1", 
              destfile = "SRR11397715.fastq.gz")


# Uncompress fasta.gz & fastq.gz files
#Rsamtools::gunzip("SRR11397715.fastq.gz", remove = FALSE)
system("time gunzip -k SRR11397715.fastq.gz")
system("time gunzip Saccharomyces_cerevisiae_genome.fa.gz")

# make some tidy 
dir.create("mv *.fa data/")
system("mv *.fastq.gz *.fastq data")




# ## ---fasta------------------------------------------------------------

# # Read fasta : Saccharomyces genome
# fa_sample <- readFasta(dirPath = "data", pattern = "genome")
# 
# # print sample
# fa_sample
# 
# # methods accessors 
methods(class = "ShortRead")
id(fa_sample)
length(fa_sample)
width(fa_sample)
sread(fa_sample)
# 

# Extract sequence from DNAStringSet object
#toString(sread(fa_sample)[1])
#as.character(sread(fa_sample)[1])

# # Write a ShortRead object 
# system("mkdir results")
# writeFasta(object = fa_sample, file = "results/my_chromosomes.fasta")




## ---fastq---------------------------------------------------------------------

# read fastq -------------------------------------------------------------------
fqsample <- readFastq(dirPath = "data", pattern = "fastq")
fqsample <- readFastq(dirPath = "data/SRR11397715.fastq")

# print fqsample
fqsample



## accesors : fastq ------------------------------------------------------------
methods(class = "ShortReadQ")
length(fqsample)
sread(fqsample)[1:20]
width(fqsample)[1:20]
sum(width(fqsample)) / 1e6 # comparar con slide "entrez" : 96.3M


hist(width(fqsample))
abline(v = 150, col = "red", lwd = 2)

fqsample150 <- fqsample[width(fqsample) >= 150]
hist(width(fqsample150))

length(fqsample) / 1e3
length(fqsample150) / 1e3



# # Write a ShortRead object 
writeFasta(object = fqsample150, file = "results/reads150.fastq")





## subset from a fastq file -----------------------------------------------------
set.seed(2022)
fqsampler <- FastqSampler("data/SRR11397715.fastq", 300000)

# extract the sample from the stored file 
fqsubset <-  yield(fqsampler)

#Warning message:
#  cerrando la conenexion 3 (reads/SRR11397715.fastq) que no esta siendo utilizada  
close(fqsampler)

length(fqsample)
length(fqsubset)

width(fqsubset)[1:20]
sread(fqsubset)[1:20]
sum(width(fqsubset))/1e6 # N Mb representados en ese subset


# otra forma : filtro 500 al azar del set de reads con 150 bp. 
length(fqsample150)
idx = sample(1:length(fqsample150), 500, replace = FALSE)
fq_sample150_500 <- fqsample150[idx]

length(fq_sample150_500)
sread(fq_sample150_500)





# ## convertQuality---------------------------------------------------------------
qreads <- quality(fqsubset) # quality encoded into ASCII characters
qreads
encoding(qreads)

dim(as(qreads , "matrix"))
as(qreads , "matrix")[1:6, 1:15]
sort(encoding(qreads)[c("A", "E", "/", "<")])


# PhredQuality instance
pq <- PhredQuality(qreads)

# transform the encoding into scores
qs <- as(pq, "IntegerList")
qs[[1]]
qs[1]
# 
# 
# # Remember : a score of 30 is considered a good quality as it means 
# # the accuracy of base call is 99.9%
# table(qs[[1]] >= 30)
# 
# # Quaility assesment & summary 
# qaSummary <- qa(fqsubset, lane = 1) 
# 
# class(qaSummary)
# methods(class = "ShortReadQQA")
# names(qaSummary)
# 
# # qa elements are accessed with qa[["name"]]
# qaSummary[["readCounts"]]
# qaSummary[["baseCalls"]]
# qaSummary[["adapterContamination"]]
# 
# 
# # get html report
# browseURL(report(qaSummary))



## Duplicates   ---------------------------------------------------------------
table(srduplicated(fqsample)) # counts of the duplicated reads 

# Cleaning reads from duplicates
uniqueReads <- fqsample[srduplicated(fqsample) == FALSE]
table(srduplicated(uniqueReads))



## Write a ShortRead object ----------------------------------------------------
# system("mkdir results")
# writeFasta(object = fq_sample150_500, file = "results/my_fasta.fasta")
# writeFastq(object = fq_sample150_500, file = "results/fqsample.fastq.gz")







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
# dir.create("data")   #(To create a directory/folder in the indicated path)
# dir("data")          #(To see files and folders in the indicated path)



# Download data ----------------------------------------------------------------
file.exists("NGS_MB_SRR9336468_1.fastq.gz")
file.exists("NGS_MB_SRR9336468_2.fastq.gz")
file.exists("Saccharomyces_cerevisiae_genome.gff3.gz")
file.exists("Saccharomyces_cerevisiae_genome.fa.gz")


download.file("https://www.dropbox.com/s/v06um7vt9ojdf42/NGS_MB_SRR9336468_1.fastq.gz?dl=1", 
              destfile = "data/1MB_SRR9336468_1.fastq.gz")

download.file("https://www.dropbox.com/s/kgxfth8ra675ccu/NGS_MB_SRR9336468_2.fastq.gz?dl=1", 
              destfile = "data/1MB_SRR9336468_2.fastq.gz")

download.file("https://www.dropbox.com/s/qtaret1hrbvw2xb/Saccharomyces_cerevisiae_genome.gff3.gz?dl=1", 
              destfile = "data/Saccharomyces_cerevisiae_genome.gff3.gz")

download.file("https://www.dropbox.com/s/4ft480eky7kghzw/Saccharomyces_cerevisiae_genome.fa.gz?dl=1", 
              destfile = "data/Saccharomyces_cerevisiae_genome.fa.gz")



# # descarga original CSIC : 
# 
# download.file(url = "https://bioinfogp.cnb.csic.es/courses/quedateencasa/23_03_2020/1M_SRR9336468_1.fastq.gz", 
#                destfile = "alignment/data/1MB_SRR9336468_1.fastq.gz")
# 
# download.file(url = "https://bioinfogp.cnb.csic.es/courses/quedateencasa/23_03_2020/1M_SRR9336468_2.fastq.gz", 
#               destfile = "alignment/data/1MB_SRR9336468_2.fastq.gz")
# 
# download.file(url = "https://bioinfogp.cnb.csic.es/courses/quedateencasa/23_03_2020/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz",
#               destfile = "alignment/data/Saccharomyces_cerevisiae_genome.fa.gz")
# 
# download.file(url = "https://bioinfogp.cnb.csic.es/courses/quedateencasa/23_03_2020/Saccharomyces_cerevisiae.R64-1-1.99.gff3.gz", 
#               destfile = "alignment/data/Saccharomyces_cerevisiae_genome.gff3.gz")



system("ls -lh")
dir()
list.files()


# Uncompress genome fastq gz files with "gunzip"
#Rsamtools::gunzip("Saccharomyces_cerevisiae_genome.gff3.gz", remove = FALSE)
system("time gunzip -k data/Saccharomyces_cerevisiae_genome.fa.gz")
system("time gunzip -k data/Saccharomyces_cerevisiae_genome.gff3.gz")

# Uncompress reads fastq gz files with "gunzip"
system("time gunzip -k data/NGS*.fastq.gz")

# Uncompress all fastq gz files with "gunzip"
system("time gunzip -k data/*.fastq.gz")


# Make some tidy
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

# Look at one of those two .fastq files
system("wc -l data/fastq/1MB_SRR9336468_1.fastq")
system("less -S fastq/1MB_SRR9336468_1.fastq")

# Look at the two reads from the same cDNA fragment 
system("head -n 4 data/fastq/1MB_SRR9336468_1.fastq")
system("head -n 4 data/fastq/1MB_SRR9336468_2.fastq")


# Notice : 
# name 1st sequence in the 1st file = name 1st seq. in the 2nd file, except 
# 1/2 and 1/2 => it is the 1st pair and the 2nd pair of the same cDNA fragment. 


# remember : the reads can contain junctions.
# muy buena explicacion Michael Love min 5 "First look at a FASTQ file" video
# ahi explica q el mismo nombre de seq aparece en los 2 ficheros

# preguntas 2ejercicio clase : "FASTQ assessment" EDx



## BLAST genome ------------------------------------------------------------
# blast genome : saccharomices
# https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&USER_FORMAT_DEFAULTS=on&SET_SAVED_SEARCH=true&PAGE=MegaBlast&PROGRAM=blastn&GAPCOSTS=0%200&MATCH_SCORES=1,-2&BLAST_SPEC=Assembly&DATABASE=genomic/559292/GCF_000146045.2&BLAST_PROGRAMS=megaBlast&MAX_NUM_SEQ=100&SHORT_QUERY_ADJUST=on&EXPECT=0.05&WORD_SIZE=28&REPEATS=4932&TEMPLATE_TYPE=0&TEMPLATE_LENGTH=0&FILTER=L&FILTER=R&FILTER=m&EQ_MENU=Enter%20organism%20name%20or%20id--completions%20will%20be%20suggested&PROG_DEFAULTS=on&SHOW_OVERVIEW=on&SHOW_LINKOUT=on&ALIGNMENT_VIEW=Pairwise&MASK_CHAR=2&MASK_COLOR=1&GET_SEQUENCE=on&NUM_OVERVIEW=100&DESCRIPTIONS=100&ALIGNMENTS=100&FORMAT_OBJECT=Alignment&FORMAT_TYPE=HTML
# blast result analysis : hacer primero @SRR9336468.1 1/2

# @SRR9336468.1 1/2 : https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=6ZKY810D013
# read mapea sitio único XV : genome data viewer
# mapea : 780 - 928
# orientación : forward
# 1a base de esa read mapea en la coord. 780 del chr XVI

# @SRR9336468.1 1/1 : https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=6ZKPUBMW013
# read mapea sitio único XV : genome data viewer
# mapea : 803 - 951
# orientación : reverse
# 1a base de esa read mapea en la coord. 951 del chr XVI








## Genome alignment ------------------------------------------------------------
#bowtie2 paper : https://www.nature.com/articles/nmeth.1923
#bowtie2 software homepage : http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
#BiocManager::install("Rbowtie2")
#library(Rbowtie2)

# STEP 1: Create genome index for bowtie2 with "bowtie2_build":
start_time = Sys.time()
bowtie2_build(references = "data/chromosomes/Saccharomyces_cerevisiae_genome.fa",
              bt2Index = "data/chromosomes/Scerevisiae_genome",
              overwrite = TRUE)
end_time = Sys.time()
end_time - start_time

# Time difference of 18.45692 secs


# bowtie2 no me genera el archivo de indice .fai que deberia haber en chromosomes/
# pero encuentro como hacerlo con Rsamtools
# samtools faidx fastafile   # samtools en linux
indexFa("data/chromosomes/Saccharomyces_cerevisiae_genome.fa") # crea el .fai


# STEP 2: Align FASTQ files against indexed genome with "bowtie2":
start_time = Sys.time()
bowtie2(bt2Index = "data/chromosomes/Scerevisiae_genome",
        samOutput = "SRR9336468.sam",
        seq1 = "data/fastq/1MB_SRR9336468_1.fastq",
        seq2 = "data/fastq/1MB_SRR9336468_2.fastq",
        overwrite = TRUE,
        "--threads=3")
end_time = Sys.time()
end_time - start_time

# Time difference of 4.716987 mins        ##"--threads=2" en despacho
# Time difference of 2.249662 mins        ##"--threads=3" en despacho
# Time difference of 1.685986 mins        ##"--threads=3" en casa


# STEP 3 : Convert SAM files into BAM (and indexes .bai) with "asBam"
# converting SAM file into BAM file (more compact and indexed)
#BiocManager::install("Rsamtools")
#library(Rsamtools)

start_time = Sys.time()
asBam("SRR9336468.sam")
end_time = Sys.time()
end_time - start_time
# Time difference of 1.345935 mins 


# make some tidy
rm(end_time, start_time)
# una vez se tiene el BAM se puede borrar el SAM porque son los mas prescindibles y los que mas ocupan
file.remove("SRR9336468.sam")                  


# Nota para mi (Curso con Esteban): GenomicAlignments package. 
# READING FROM A BAM FILE 
# part 1 : https://www.youtube.com/watch?v=2KqBSbkfhRo
# part 2 : https://www.youtube.com/watch?v=3PK_jx44QTs&t=140s
# library(GenomicAlignments)
# gal <- readGAlignments("SRR9336468.bam")
# gal
# length(gal)
# table(strand(gal))




## IGV  ----
#ver video Michael Love 

# IGV software homepage : https://www.broadinstitute.org/igv/
# https://igv.org/
  


# Load genome : Genome / Local File / chromosomes / cargar a la vez el .fa + .fai
# Load genes : Tracks / Local File (gff3)
#View options : 
  #Set track height / 100
  #Expanded / Collapse
# Load reads : Tracks / Load File / select bam + index (bai)


  # select chr1
#  - note you have reads where you have genes. 
#  - note some genes have more reads than others = in the experimental conditions , 
#    some genes are highly expressed + some genes are not expressed. 
#   - explicar mistmach 
  #        - small variation related to the reference genome 
  #      - (n. mismatch are allowed during alignment step)
  #      - if you find a consistent mistmatch (several times on the same column) : 
  #          - the reference genome has error
  #          - the cell contains a mutation (SNP)
  #         - 57,300 - 58,000 : YAL045C

# color by / read strand 
# especificidad de hebra : color by / first-of-pair strand
#    5'-3' : azul
#    3'-5' : naranja
  
  
  
# preguntas 2ejercicio clase : "IGV assessment" 





## Filtering reads by mapping quality
# ------------------------------------------------------------------------------
bamfile <- BamFile("SRR9336468.bam")
bamfile
seqinfo(bamfile)

# set the filter condition
# build the filtered : .bam + .bai files 
param = ScanBamParam(mapqFilter = 35) # min=35 * 
dest = filterBam(file = bamfile, destination = "results/mapq.bam", param = param) 

#Tb puedo controlar la calidad del alineamiento despues en Deseq2 con la
#function 'featureCounts' y la variable minMQS

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


# FIN CLASE






# Nota para mi : No sé filtrar duplicados PCR/opticos porque si 
# defino el filtro expresamente, el numero de lecturas es el mismo que si no 
# lo hago. Y si pongo TRUE sigue saliendo el mismo n. lecturas 

param2 = ScanBamParam(mapqFilter = 35, 
                      flag = scanBamFlag(isDuplicate = FALSE)) 

dest2 = filterBam(file = bamfile, 
                  destination = "results/mapqNOdup.bam", 
                  param = param2) 


#dest2 <- BamFile("results/mapqNOdup.bam")

countBam(bamfile)
countBam(dest)
countBam(dest2)




# IGV only with reads >= 35 mapping quality score
# comparar diferencias : chr1 _ 78kb _ FUN12







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


# Paper : 
# "Physiological responses of Saccharomyces cerevisiae to industrially
# relevant conditions: Slow growth, low pH, and high CO2 levels"
# Hakkaart X, Liu Y, Hulst M, El Masoudi A, Peuscher E, Pronk J, van Gulik W, Daran-Lapujade P.
# Biotechnol Bioeng. 2020 Mar;117(3):721-735. doi: 10.1002/bit.27210. Epub 2020 Jan 22.

# https://www.ncbi.nlm.nih.gov/pubmed/31654410


# sample table :
# https://github.com/jdieramon/master_biotec_NGS/blob/main/sample_table.md

# SampleID      Description               Replicate
# SRR9336468    Chemostat pH5 0.04% CO2   R1
# SRR9336469    Chemostat pH5 0.04% CO2   R2
# SRR9336470    Chemostat pH5 0.04% CO2   R3
# SRR9336471    Chemostat pH5 50% CO2     R1
# SRR9336472    Chemostat pH5 50% CO2     R2
# SRR9336473    Chemostat pH5 50% CO2     R3
# SRR9336474    Chemostat pH3 0.04% CO2   R1
# SRR9336475    Chemostat pH3 0.04% CO2   R2
# SRR9336476    Chemostat pH3 0.04% CO2   R3

# install.packages("BiocManager") # <- Bioconductor manager

# Install R packages (from Bioconductor):

# BiocManager::install("Rbowtie2")  # <- align FASTQ files
# BiocManager::install("Rsamtools") # <- manipulate SAM files
# BiocManager::install("Rsubread")  # <- quantify alignments
# BiocManager::install("DESeq2")    # <- differential gene expression

# load all libraries:

#library("R.utils")
library("Rbowtie2")
library("Rsamtools")
library("Rsubread")
library("DESeq2")



## STEP 1: Download gz files 
# -----------------------------------------------------------------------------


# Download genome data
download.file("https://www.dropbox.com/s/4ft480eky7kghzw/Saccharomyces_cerevisiae_genome.fa.gz?dl=1", 
              destfile = "Saccharomyces_cerevisiae_genome.fa.gz")

download.file("https://www.dropbox.com/s/qtaret1hrbvw2xb/Saccharomyces_cerevisiae_genome.gff3.gz?dl=1", 
              destfile = "Saccharomyces_cerevisiae_genome.gff3.gz")



# Download FASTQ files (HUGE!!!)
# https://www.ebi.ac.uk/ena/data/view/SRR9336468
# https://www.ebi.ac.uk/ena/data/view/SRR9336469
# https://www.ebi.ac.uk/ena/data/view/SRR9336470
# https://www.ebi.ac.uk/ena/data/view/SRR9336471
# https://www.ebi.ac.uk/ena/data/view/SRR9336472
# https://www.ebi.ac.uk/ena/data/view/SRR9336473
# https://www.ebi.ac.uk/ena/data/view/SRR9336474
# https://www.ebi.ac.uk/ena/data/view/SRR9336475
# https://www.ebi.ac.uk/ena/data/view/SRR9336476

# fastq.gz = 1.6 GB
# fastq = 7.4 GB


#forma 1 (fichero original ; me da error)
download.file("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR933/008/SRR9336468/SRR9336468_1.fastq.gz", 
              destfile = "SRR9336468_1.fastq.gz")

#forma2 (fichero original)
url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR933/008/SRR9336468/SRR9336468_1.fastq.gz"
curl::curl_download(url, destfile = "SRR9336468_1.fastq.gz")



# library(ShortRead)
# 
# # read fastq ------------------------------------------------------------------
# start_time = Sys.time()
# faq_sample <- readFastq(dirPath = "SRR9336468_1.fastq")
# end_time = Sys.time()
# end_time - start_time
# #Time difference of 4.461314 mins
# 
# faq_sample
# 
# ## accesorFastq ---------------------------------------------------------------
# methods(class = "ShortReadQ")
# id(faq_sample)
# length(faq_sample)
# sread(faq_sample)
# width(faq_sample)
# table(width(faq_sample))

# rm(faq_sample)

# dropbox (1M reads) - en clase descargar este. 
download.file("https://www.dropbox.com/s/v06um7vt9ojdf42/NGS_MB_SRR9336468_1.fastq.gz?dl=1", 
              destfile = "1M_SRR9336468_1.fastq.gz")

download.file("https://www.dropbox.com/s/kgxfth8ra675ccu/NGS_MB_SRR9336468_2.fastq.gz?dl=1", 
              destfile = "1M_SRR9336468_2.fastq.gz")


# no es necesario que descarguen el resto de fastq.gz files porque les dare 
# despues los bam files 




## STEP 2: Uncompress all gz files with "gunzip"
# -----------------------------------------------------------------------------

# Uncompress genome 
#Rsamtools::gunzip("Saccharomyces_cerevisiae_genome.gff3.gz", remove = FALSE)
system("time gunzip -k Saccharomyces_cerevisiae_genome.fa.gz")
system("time gunzip -k Saccharomyces_cerevisiae_genome.gff3.gz")

dir.create("genes")
system("mv Saccharomyces_cerevisiae_genome.* genes")



# Uncompress reads
dir.create("fastq")
system("time gunzip -k *.fastq.gz")
system("mv *.fastq *.gz fastq/")





## STEP 3: Create genome index for bowtie2 with "bowtie2_build":
# ----------------------------------------------------------------------------- 
bowtie2_build(references = "genes/Saccharomyces_cerevisiae_genome.fa", 
              bt2Index = "genes/Sc_genome")

# bowtie2 no me genera el archivo de indice .fai
# pero encuentro como hacerlo con Rsamtools
# samtools faidx some.fasta   # samtools en linux
indexFa("genes/Saccharomyces_cerevisiae_genome.fa") # crea el .fai





# STEP 4: Align FASTQ files against indexed genome with "bowtie2"
# ----------------------------------------------------------------------------- 

dir.create("bam")

start_time = Sys.time()
bowtie2(bt2Index = "genes/Sc_genome", 
        samOutput = "bam/1M68_pH5_0.04CO2_R1.sam", 
        seq1 = "fastq/1M_SRR9336468_1.fastq", 
        seq2 = "fastq/1M_SRR9336468_2.fastq", "--threads=3")       
end_time = Sys.time()
end_time - start_time
# Time difference of 2.249662 mins        ##"--threads=3" en despacho


# bowtie2(bt2Index = "genes/Sc_genome", samOutput = "bam/1M69_pH5_0.04CO2_R2.sam", seq1 = "fastq/1M_SRR9336469_1.fastq", seq2 = "fastq/1M_SRR9336469_2.fastq", "--threads=3")
# bowtie2(bt2Index = "genes/Sc_genome", samOutput = "bam/1M70_pH5_0.04CO2_R3.sam", seq1 = "fastq/1M_SRR9336470_1.fastq", seq2 = "fastq/1M_SRR9336470_2.fastq", "--threads=3")
# 
# bowtie2(bt2Index = "genes/Sc_genome", samOutput = "bam/1M71_pH5_5CO2_R1.sam", seq1 = "fastq/1M_SRR9336471_1.fastq", seq2 = "fastq/1M_SRR9336471_2.fastq", "--threads=3")
# bowtie2(bt2Index = "genes/Sc_genome", samOutput = "bam/1M72_pH5_5CO2_R2.sam", seq1 = "fastq/1M_SRR9336472_1.fastq", seq2 = "fastq/1M_SRR9336472_2.fastq", "--threads=3")
# bowtie2(bt2Index = "genes/Sc_genome", samOutput = "bam/1M73_pH5_5CO2_R3.sam", seq1 = "fastq/1M_SRR9336473_1.fastq", seq2 = "fastq/1M_SRR9336473_2.fastq", "--threads=3")
# 
# bowtie2(bt2Index = "genes/Sc_genome", samOutput = "bam/1M74_pH3_0.04CO2_R1.sam", seq1 = "fastq/1M_SRR9336474_1.fastq", seq2 = "fastq/1M_SRR9336474_2.fastq", "--threads=3")
# bowtie2(bt2Index = "genes/Sc_genome", samOutput = "bam/1M75_pH3_0.04CO2_R2.sam", seq1 = "fastq/1M_SRR9336475_1.fastq", seq2 = "fastq/1M_SRR9336475_2.fastq", "--threads=3")
# bowtie2(bt2Index = "genes/Sc_genome", samOutput = "bam/1M76_pH3_0.04CO2_R3.sam", seq1 = "fastq/1M_SRR9336476_1.fastq", seq2 = "fastq/1M_SRR9336476_2.fastq", "--threads=3")

 


# # STEP 5: Convert SAM files into BAM (and indexes .bai) with "asBam"
# ----------------------------------------------------------------------------- 
asBam(file = "bam/1M68_pH5_0.04CO2_R1.sam", destination = "bam/basal_1")
# asBam("bam/1M69_pH5_0.04CO2_R2.sam", "bam/basal_2")
# asBam("bam/1M70_pH5_0.04CO2_R3.sam", "bam/basal_3")
# 
# asBam("bam/1M71_pH5_5CO2_R1.sam", "bam/indust1_1")
# asBam("bam/1M72_pH5_5CO2_R2.sam", "bam/indust1_2")
# asBam("bam/1M73_pH5_5CO2_R3.sam", "bam/indust1_3")
# 
# asBam("bam/1M74_pH3_0.04CO2_R1.sam", "bam/indust2_1")
# asBam("bam/1M75_pH3_0.04CO2_R2.sam", "bam/indust2_2")
# asBam("bam/1M76_pH3_0.04CO2_R3.sam", "bam/indust2_3")


## una vez se tiene el BAM se puede borrar el SAM porque son los mas prescindibles y los que mas ocupan

junk <- list.files(path = "bam", pattern = "*.sam", full.names = TRUE)
file.remove(junk)



# download the rest of *.bam files
start_time = Sys.time()
download.file("https://www.dropbox.com/s/pr3vqfj9r2c35xg/bam_files.zip?dl=1", 
              destfile = "bam/bam_files.zip")
end_time = Sys.time()
end_time - start_time
# Time difference of 1.000254 mins (casa)
# Time difference of 2.195933 mins (despacho) 



# uncompress bam files
unzip(zipfile = "bam/bam_files.zip", exdir = "bam")
# si da error, es porque la descarga del fichero .zip no ha funcionado bien. 
# el tamaño del archivo  correcto es 871MB
# si se repite y sigue sin funcionar, que copien la dirección en el navegador internet
# se descarga el archivo y lo mueven a la carpeta 

file.remove("bam/bam_files.zip")




## STEP 7: store all paths of bam files into "bamFiles" R object (use "list.files"):
# ----------------------------------------------------------------------------- 
list.files(path = "bam", pattern = "*.bam$", full.names = TRUE)
bamFiles <- list.files(path="bam", pattern="*.bam$", full.names=TRUE)







# STEP 5: Create metadata file for the Saccharomyces count matrix
# ----------------------------------------------------------------------------- 
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






# # STEP 8: quantify "genes" with "featureCounts" and store results in R object "data"
# -----------------------------------------------------------------------------  
# # [OPTIONAL] To check the content of gff3 file:
#file.show("genes/Saccharomyces_cerevisiae_genome.gff3")

start_time = Sys.time()
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
  nthreads = 4          # by default = 1.
)
end_time = Sys.time()
end_time - start_time

# Time difference of 11.27765 secs (threads = 3)
# Time difference of 44.48049 secs (threads = 3)
 
# class(data)
# data$counts[1:10, 1:3]
# data$counts = raw_counts

dir.create("results")
write.table(data$counts, file = "results/raw_counts" , sep = "\t")

# .tsv means "tab separated values"
raw_counts = read.table("results/raw_counts")


# trabajar con la clase los siguientes conceptos : 
# mapping gen YAL063C in the last 4 samples 
raw_counts["gene:YAL063C", 6:9]

# does the mapping match what I see in IGV ?


# open IGV : https://igv.org/
# load Saccharomices genome
# search : YAL063C (v. desktop no permite busqueda ; app sí lo encuentra y 
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
abline(0,1, col = "red")


plot(log.norm.counts[,c(1,7)], cex=.1)
abline(0,1, col = "red")


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
samples
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

# recordar que estos genes 'sig_norm_counts' salen de haber hecho el contraste
# basal vs ind1. Es normal que el heatmap recoja esta diferencia. 


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

