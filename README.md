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
BiocManager::install("Rbowtie2")  # <- align FASTQ files
BiocManager::install("Rsamtools") # <- manipulate SAM files
BiocManager::install("Rsubread")  # <- quantify alignments
BiocManager::install("DESeq2")    # <- differential gene expression
```
