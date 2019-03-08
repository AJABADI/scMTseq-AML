options(warn=-1)
suppressMessages(library(data.table))
# suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressMessages(library(seqinr))
suppressMessages(library(stringr))


io$basedir <- "../data/scBSseq"
io$anno.folder <- str_c(io$basedir,"/filt")

io$in.folder <- str_c(io$basedir,"/met/raw/bismark")
io$out.folder <- str_c(io$basedir,"/met/parsed")

## valid chromosomes
opts$chr_list <- c(1:22)

# Annotations to analyse - 'all' or names of annotation in anno.folder
opts$annos <- "all"
if (opts$annos == "all")
  opts$annos <- sapply(str_split(list.files(io$anno.folder, pattern = "\\.bed$"),"\\."),"[[", 1)

## samples
samples_keep <- sapply(str_split(list.files(io$in.folder, pattern = "\\.cov$"),"\\.cov"),"[[", 1)
stopifnot(all(!duplicated(samples_keep)))

cat(sprintf("- Processing samples"))

############################
## Preprocess annotations ##
############################

anno_list <- list()
# opts$annos = opts$annos[c(1,5)]
for (i in opts$annos) {
  
  # Read annotation file
  anno.file <- sprintf("%s/%s.bed",io$anno.folder,i)
  ## read the annotation file which is tab delimited and keep 
  dat_anno <- fread(anno.file ,sep="\t", header=F, verbose=F) %>% .[,-6] %>% 
    setnames(c("chr","start","end","strand","id"))
  
  # Check that there are no weird chromosomes
  anno_list[[i]] <- dat_anno %>% .[chr%in%opts$chr_list,] %>% setkey(chr,start,end)
}


dir.create(sprintf("%s/tmp",io$out.folder), recursive=T)


for (i in 1:length(samples_keep)) {
  sample=samples_keep[i]
  ## check if already processed
  samples_processed <- list.files(sprintf("%s/tmp",io$out.folder))
  if (sum(str_detect(samples_processed,str_c(sample,"_",opts$anno))) == length(opts$anno)) {
    cat(sprintf("Sample %s already processed for all required annotations...\n",sample)) 
  } else {
    cat(sprintf("Sample %s has not been processed, annotating...\n",sample))  
    
    # Read and parse raw methylation data (cov files)
    dat_sample <- fread(paste0(io$in.folder,'/',sample, '.cov'), sep="\t", header="auto", verbose=F, showProgress=F)
    
    # dat_sample <- dat_sample[,-3] ## drop the start which is the same as end
    ## no need for these as the cov files already have start and end
    colnames(dat_sample) <- c("chr","pos","rate","c_meth", "c_unmeth")
    # PICKUP: what is a sensible averaging?
    # Add 'start' and 'end' columns to do the overlap
    dat_sample <- dat_sample[,c("start","end") := list(pos,pos)][,pos:=NULL] %>% .[,chr:=as.factor(chr)] %>% setkey(chr,start,end)
    colnames(dat_sample) <- c("chr","start","end","rate")
    ############# must set key for foverlaps - keys are the ones that are overlapped ########
    dat_sample %<>% .[,chr:=as.factor(chr)] %>% setkey(chr,start,end)
    
    # Overlap data with annotations
    
    for (anno in opts$anno) {
      fname.out <- sprintf("%s/tmp/%s_%s.tsv",io$out.folder,sample,anno)
      if (file.exists(paste0(fname.out,".gz"))) {
        cat(sprintf("Annotation for %s with %s already found, loading...\n",sample,anno))
      } else {
        cat(sprintf("Annotating %s with %s annotations...\n",sample,anno))
        
        # Overlap - genomic context
        ov <- foverlaps(dat_sample, anno_list[[anno]], nomatch=0) %>% 
          .[,"i.end":=NULL] %>% setnames("i.start","pos") ## remove the duplicate end column and rename  i.start to pos
        
        # Calculate methylation status for each region in the annotation by summarising over all sites
        # out <- ov[,c("sample","anno") := list(sample,anno)] %>% .[,.(rate=round(mean(rate)*100), weight=.N), keyby=.(sample,id,anno)]
        out <- ov[,c("sample","anno") := list(sample,anno)] %>% .[,.(rate=round(mean(rate)), weight=.N), keyby=.(sample,id,anno)]
        
        # Store and save results
        fwrite(out, fname.out, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
      }
    }
    ######################## Windows
    
    # for (anno in paste0("anno_", opts$wins)) {
    #   fname.out <- sprintf("%s/tmp/%s_%s.tsv",io$out.folder,sample,anno)
    #   if (file.exists(paste0(fname.out,".gz"))) {
    #     cat(sprintf("Annotation for %s with %s already found, loading...\n",sample,anno))
    #   } else {
    #     cat(sprintf("Annotating %s with %s annotations...\n",sample,anno))
    #     
    #     # Overlap - genomic context
    #     ov <- foverlaps(dat_sample, get(anno), nomatch=0) %>% ## overlap them an don't retun anything if there's no match
    #       .[,"i.end":=NULL] %>% setnames("i.start","pos") ## remove the duplicate end column and rename  i.start to pos
    #     
    #     # Calculate methylation status for each region in the annotation by summarising over all sites
    #     # out <- ov[,c("sample","anno") := list(sample,anno)] %>% .[,.(rate=round(mean(rate)*100), weight=.N), keyby=.(sample,id,anno)]
    #     out <- ov[,c("sample","anno") := list(sample,anno)] %>% .[,.(rate=round(mean(rate)), weight=.N), keyby=.(sample,id,anno)]
    #     
    #     # Store and save results
    #     fwrite(out, fname.out, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
    #   }
    # }

  }
}#)

