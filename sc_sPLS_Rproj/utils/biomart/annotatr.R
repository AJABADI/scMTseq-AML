library(data.table); library(magrittr);library(dplyr);library(GenomicRanges);library(rtracklayer)
library(annotatr);library(TxDb.Hsapiens.UCSC.hg38.knownGene);library(org.Hs.eg.db)
library(BSgenome);library(BSgenome.Hsapiens.UCSC.hg38)

## all available ones for hg38:
builtin_annotations() %>% .[substring(.,0,4)=="hg38"]

Basicgenes <-"hg38_basicgenes"
Enhancer <- "hg38_enhancers_fantom"
Intergenic <- "hg38_genes_intergenic"
Splicing<- c("hg38_genes_intronexonboundaries","hg38_genes_exonintronboundaries")
CpG<- "hg38_cpgs"
CGIs<- "hg38_cpg_islands"

Basicgenes_annotation <-  build_annotations(genome='hg38',annotations =Basicgenes)
Enhancer_annotation <-  build_annotations(genome='hg38',annotations =Enhancer)
Intergenic_annotation <-  build_annotations(genome='hg38',annotations =Intergenic)
Splicing_annotation <- build_annotations(genome='hg38',annotations =Splicing)
CpG_annotation <- build_annotations(genome='hg38',annotations = CpG)
CGI_annotation <- build_annotations(genome='hg38',annotations = CGIs )

## create bed files from enhancers
## TODO What about the strands, do they matter when finding enhancers for genes?
enhancers <- granges2bed(Enhancer_annotation ) ## from utils load granges2bed
enhancers %<>% .[,1:3]
enhancers[,id:=paste("enhancer",chr,start, sep = "_")]
# save(enhancers, file="/Users/alabadi/Documents/_Projects/sc_sPLS/data/scBSseq/filt/enhancers_chr_start_end.bed")

## for every gene, look -1Mb +1Mb away from start and associate it to the closest enhancers
size=20
gene_ex <- data.table(chr=1, start = seq(10,100,10), id=paste0("id",1:10))
gene_ex
gene_ex_cenetred_1Mb <- gene_ex %>% .[,`:=`(start=start-size, end=start+size, tss=start)]
gene_ex_cenetred_1Mb[start<0,start:=0] %>% setkey(chr,start,end)
enhancer_ex <- data.table(chr=1, start = seq(20,40,10), end=seq(25,45,10), enhancer_id=paste0("enhancer-",1:3))
## long ones first, keyed, fat one second
ovlap <- foverlaps(enhancer_ex, gene_ex_cenetred_1Mb) ## overlaps with many, we want the closest
head(ovlap) ## we want the enhancer start to be close to be closest to 
## for the same overlaps (chr, start, end, i.start, i.end) we want the id with the least abs(tss-i.start)
ovlap[,dist:=abs(tss-i.start)]
## sort by decreasing dist
ovlap %<>% .[base::order(dist)]
head(ovlap)
## remove duplicates
ovlap %<>% .[!base::duplicated(ovlap, by=c("enhancer_id"))]
