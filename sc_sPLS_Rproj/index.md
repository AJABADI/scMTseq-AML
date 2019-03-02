--- 
title: "scMTseq of AML cells - pilot"
author: "Al J Abadi"
date: "2019-02-11"
bibliography:
- packages.bib
- citations.bib
site: bookdown::bookdown_site
documentclass: book
biblio-style: apalike
link-citations: no
output:
  bookdown::gitbook:
    includes:
      in_header: codefold.html
description: "Joint methylation and transcriptome sequencing of AML cells"
---



# Data {-}

Data from Dr Heather Lee lab from The University of Newcastle. 


<!-- The cells are from the Kasumi-1 AML cell line. Cell identifiers are given by a letter (indicating plate row), and a number (indicating plate column), e.g. *_A7_*, *_C10_*. -->


<!-- +  **scRNA-seq Data** -->

<!-- Count matrices which are quality controlled. -->


<!-- +  **scBS-seq Data** -->

<!-- Each .cov file in which each line is a cytosine residue, and the columns are as follows: -->
<!-- <chromosome>  <position>  <strand>  <count methylated>  <count unmethylated>  <C-context> <trinucleotide context> -->
<!-- These reports are produced by Bismark. The functional element is ill defined. Sometimes we look at promoters, sometimes enhancers, sometimes unbiased - occasionally overlapping - 3kb windows to ensure minimisation of missed genomic features. -->

<!-- 40 cells amcth those from RNA-seq. -->


```r
## update as needed.
## installing the required packages for this analysis, if necessary
required.pkgs = c('mixOmics','SingleCellExperiment','scran','data.table',
                  'DESeq2', 'edgeR', 'Rtsne', 'ggplot2', 'ggrepel', 'gridExtra',
                  'grid', 'reshape2',
                  ## for summarising
                  'seqinr', 'stringr', 'doParallel', 'argparse'
                  )
```


```r
## make sure Biocmanager is installed
if (!requireNamespace('BiocManager', quietly = T)){
  paste('Trying to install BiocManager')
  install.packages('BiocManager')
}
```


```r
## package installer function - for those not already installed
package.installer = function(pkgs=required.pkgs){
  for (package in pkgs){
    if (!requireNamespace(package, quietly = T)){
  paste0('Trying to install ', package)
  BiocManager::install(package, update = F)
    }
    }
}
## run function
package.installer(required.pkgs)
```

<!-- Activate when params are finalised and defined in YAML -->
<!-- ```{r} -->
<!-- ## I/O -->
<!-- io=list() -->
<!--   ## raw RNA data -->
<!--   io$rna_file = params$rna_file ## '../data/scRNA-seq/raw/Kasumi1_scRNAseq_RawCounts.txt' -->
<!--   ## raw single bp met calls folder -->
<!--   io$met_1bp_in = params$met_1bp_in ## '../data/scBS-seq/met/1bp/raw' -->
<!--   ## processed single bp met calls file name -->
<!--   io$met_1bp_out = params$met_1bp_out ## '../data/scBS-seq/met/1bp/parsed/met_data.tsv' -->
<!--   ## whether to use parsed data or (re-)parse -->
<!--   io$use.parsed= params$use.parsed## TRUE -->
<!--   ## 3kb windows data folder -->
<!--   io$met_3kb = params$met_3kb ## '../data/scBS-seq/met/3kbp/raw/methylVariance-w3K-s1K5.csv' -->
<!--   ## .bed files for summarising calls on genomic features -->
<!--   io$bedfiles = params$bedfiles ## '../data/scBS-seq/filt' -->
<!--   io$RData = params$RData ## save rdata -->
<!-- ``` -->



