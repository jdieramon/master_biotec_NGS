# Análisis genómicos y transcriptómicos con plataforma NGS




## Enlaces R & Bioconductor

* [Central R Archive Network (CRAN)](http://cran.rstudio.com/)
* [RStudio](http://www.rstudio.com/)
* [Bioconductor](http://bioconductor.org/install)

Una vez instalado R y el paquete `BiocManager` el siguiente codigo instalará Bioconductor:

```
BiocManager::install()
```

Para saber la version de Bioconductor que estamos usando y si los paquetes estan actualizados :  

```
BiocManager::version()
BiocManager::valid()
```

## Instalación paquetes Bioconductor
```
BiocManager::install("ShortRead")  # <- manipulate FASTQ files
BiocManager::install("Rbowtie2")  # <- align FASTQ files
BiocManager::install("Rsamtools") # <- manipulate SAM files
BiocManager::install("Rsubread")  # <- quantify alignments
BiocManager::install("DESeq2")    # <- differential gene expression
```

## Material 
  * Análisis expresion Diferencial RNA-Seq : [tabla muestras](sample_table.md)
  * European Nucleotide Archive  
    [Run: SRR9336468](https://www.ebi.ac.uk/ena/browser/view/SRR9336468)  
    [Run: SRR9336469](https://www.ebi.ac.uk/ena/browser/view/SRR9336469)  
    [Run: SRR9336470](https://www.ebi.ac.uk/ena/browser/view/SRR9336470)   
    [Run: SRR9336471](https://www.ebi.ac.uk/ena/browser/view/SRR9336471)   
    [Run: SRR9336472](https://www.ebi.ac.uk/ena/browser/view/SRR9336472)   
    [Run: SRR9336473](https://www.ebi.ac.uk/ena/browser/view/SRR9336473)   
    [Run: SRR9336474](https://www.ebi.ac.uk/ena/browser/view/SRR9336474)   
    [Run: SRR9336475](https://www.ebi.ac.uk/ena/browser/view/SRR9336475)   
    [Run: SRR9336476](https://www.ebi.ac.uk/ena/browser/view/SRR9336476) 
    


## Descarga ficheros secuencia genómica         
```
download.file("https://www.dropbox.com/s/4ft480eky7kghzw/Saccharomyces_cerevisiae_genome.fa.gz?dl=1", 
              destfile = "Saccharomyces_cerevisiae_genome.fa.gz")
````

## Descarga ficheros anotación genómica         
```
download.file("https://www.dropbox.com/s/qtaret1hrbvw2xb/Saccharomyces_cerevisiae_genome.gff3.gz?dl=1", 
              destfile = "Saccharomyces_cerevisiae_genome.gff3.gz")              
```

## Descarga ficheros lecturas
```
download.file("https://www.dropbox.com/s/v06um7vt9ojdf42/NGS_MB_SRR9336468_1.fastq.gz?dl=1", 
              destfile = "1M_SRR9336468_1.fastq.gz")

download.file("https://www.dropbox.com/s/kgxfth8ra675ccu/NGS_MB_SRR9336468_2.fastq.gz?dl=1", 
              destfile = "1M_SRR9336468_2.fastq.gz")
````

## Mapping Reads
   
 * STEP 1: Create genome index for bowtie2 with `bowtie2_build` 
 * STEP 2: Align FASTQ files against indexed genome with `bowtie2` 
 * STEP 3: Convert SAM files into BAM (and indexes .bai) with `asBam` 

[code](mappingReads.R)
