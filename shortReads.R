
###########################################################
#                                                         #
#  Sequence reads assessment
#                                                         #
#                                                         #
#                                                         #
###########################################################

# Analisis de genomas y transcriptomas con plataforma NGS
# Master Biotecnologia
# APR 2023

# Jose V. Die, Dep Genetics, UCO
# jose.die@uco.es


# Install R packages (from CRAN):
install.packages("BiocManager").   # <- Bioconductor manager
#install.packages("R.utils")       # <- uncompress *.gz files

# Install R packages (from Bioconductor):
BiocManager::install("ShortRead")  # <- manipulate *.fastq files

# Dependencies
library(ShortRead)  


# Organize the Project
dir.create("data")
dir.create("results")


# Download files 
#https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_000146045.2/
download.file("https://www.dropbox.com/s/4ft480eky7kghzw/Saccharomyces_cerevisiae_genome.fa.gz?dl=1", 
              destfile = "data/Saccharomyces_cerevisiae_genome.fa.gz")

#https://www.ncbi.nlm.nih.gov/sra/?term=SRR11397715
download.file(url = "https://www.dropbox.com/s/svch71j19zfhntm/SRR11397715.fastq.gz?dl=1", 
              destfile = "data/SRR11397715.fastq.gz")


# Uncompress fasta.gz & fastq.gz files
#Rsamtools::gunzip("SRR11397715.fastq.gz", remove = FALSE)
system("time gunzip -k data/SRR11397715.fastq.gz")
system("time gunzip -k data/Saccharomyces_cerevisiae_genome.fa.gz")



# ## ---fasta------------------------------------------------------------------

# Read fasta : Saccharomyces genome
fa_sample <- readFasta(dirPath = "data/Saccharomyces_cerevisiae_genome.fa")

# Print sample
fa_sample

# Accesors ShortRead
methods(class = "ShortRead")
id(fa_sample)
length(fa_sample)
width(fa_sample)
sread(fa_sample)

# Extract sequence from DNAStringSet object
toString(sread(fa_sample)[1])
as.character(sread(fa_sample)[1])

# Write a ShortRead object 
writeFasta(object = fa_sample, file = "results/Saccchromosomes.fasta", compress = TRUE)



## ---fastq---------------------------------------------------------------------

# Read fastq 
fqsample <- readFastq(dirPath = "data", pattern = "fastq")
fqsample <- readFastq(dirPath = "data/SRR11397715.fastq")

# Print fqsample
fqsample


# Accesors : fastq 
methods(class = "ShortReadQ")
length(fqsample)
sread(fqsample)[1:20]
width(fqsample)[1:20]
sum(width(fqsample)) / 1e6 # comparar con slide "entrez" : 96.3M


# Plotting
hist(width(fqsample))
abline(v = 150, col = "red", lwd = 2)

fqsample150 <- fqsample[width(fqsample) >= 150]
hist(width(fqsample150))

length(fqsample) / 1e3
length(fqsample150) / 1e3


# Write a ShortRead object 
writeFasta(object = fqsample150, file = "results/reads150.fastq")




## subset from a fastq file -----------------------------------------------------
set.seed(2023)
fqsampler <- FastqSampler("data/SRR11397715.fastq", 300000)

# extract the sample from the stored file 
fqsubset <-  yield(fqsampler)

#Warning message:
#  cerrando la conexion (reads/SRR11397715.fastq) que no esta siendo utilizada  
close(fqsampler)

length(fqsample)
length(fqsubset)

width(fqsubset)[1:20]
sread(fqsubset)[1:20]
sum(width(fqsubset))/1e6 # N Mb representados en ese subset


# alternative: filter out 500 random sequences fron the 150bp dataset. 
length(fqsample150)
sample(1:10, 3, replace = FALSE)
set.seed(2023)
idx = sample(1:length(fqsample150), 500, replace = FALSE)
fq_sample150_500 <- fqsample150[idx]

length(fq_sample150_500)
sread(fq_sample150_500)



## convertQuality---------------------------------------------------------------
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
 
 
# Remember : a score of 30 is considered a good quality as it means
# the accuracy of base call is 99.9%
table(qs[[1]] >= 30)

# Quaility assesment & summary
qaSummary <- qa(fqsubset, lane = 1)

class(qaSummary)
methods(class = "ShortReadQQA")
names(qaSummary)

# qa elements are accessed with qa[["name"]]
qaSummary[["readCounts"]]
qaSummary[["baseCalls"]]
qaSummary[["adapterContamination"]]


# get html report
browseURL(report(qaSummary))



## Duplicates   ---------------------------------------------------------------
table(srduplicated(fqsample)) # counts of the duplicated reads 

# Cleaning reads from duplicates
uniqueReads <- fqsample[srduplicated(fqsample) == FALSE]
table(srduplicated(uniqueReads))



## Write a ShortRead object ----------------------------------------------------
# system("mkdir results")
# writeFasta(object = fq_sample150_500, file = "results/my_fasta.fasta")
# writeFastq(object = fq_sample150_500, file = "results/fqsample.fastq.gz",compress=TRUE)
