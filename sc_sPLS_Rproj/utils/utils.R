############################ global plot helper
anno_plot_vec <- function(dt=pearson_value){
  contexts <-  c(genebody="Genebody", prom="Promoter",prom_cgi="Promoter-CGI", prom_noncgi="Promoter-nonCGI", intergenic="Intergenic", CGI="CG Islands")
  colz <- gg_color_hue(length(contexts)) %>% set_names(names(contexts))
  
  contexts <- contexts %>% .[names(.) %in% as.character(unique(dt$anno))]
  colz <- colz %>% .[names(.) %in% as.character(unique(dt$anno))]
  
  anno_vec <-  colz %>% structure(., levels=unname(contexts))
  return(anno_vec)
}
############################ gg box plot by annotation - uses a df_anno data.frame
gg_box <- function(dt,aesX,aesY, global=NULL,...){
  
  anno_vec <- anno_plot_vec(dt)
  
  dt$anno %<>% factor(., levels = names(anno_vec), ordered = TRUE)
  dt %<>% .[order(dt$anno, decreasing = FALSE),]
  
  p <- ggplot(dt) 
  if(!is.null(global)){
    # p <- p + geom_rect(aes(xmin=min(dt[,aesX, with=F]), xmax=max(dt[,aesX, with=F]), ymin=global[1], ymax=global[6]), fill="grey80")
    # p <- p + geom_rect(aes(xmin=min(dt[,aesX, with=F]), xmax=max(dt[,aesX, with=F]), ymin=global[2], ymax=global[5]), fill="grey60")
    p <- p + geom_rect(aes(xmin=0, xmax=length(anno_vec)+0.5, ymin=global[1], ymax=global[6]), fill="grey80")
    p <- p + geom_rect(aes(xmin=0, xmax=length(anno_vec)+0.5, ymin=global[2], ymax=global[5]), fill="grey60")
  }
  
  p <- p + theme_bw() + 
    geom_boxplot(aes_string(x=aesX, y=aesY,  col=aesX),alpha=0.7, show.legend = FALSE, width=0.4) +labs(...) +
    geom_violin(aes_string(x=aesX, y=aesY, fill=aesX, col=aesX),alpha=0.7, show.legend = FALSE) +
    scale_fill_manual(values = anno_vec, label=levels(anno_vec), guide_legend(title="Context")) +
    scale_color_manual(values = anno_vec, label=levels(anno_vec), guide_legend(title="Context")) +
    scale_x_discrete(labels=levels(anno_vec)) + theme(axis.text.x = element_text(face="bold")) + coord_flip()
  # scale_x_discrete(position = "top")
  
  return(p+theme(axis.text.y = element_text(colour="black")))
}

# ############################ global plot helper
# anno_df <- function(annos = NULL,
#                           contexts = NULL){
#   if(is.null(annos) | is.null(contexts)){
#     contexts <- c("Genebody", "Intergenic", "Promoter","Promoter-CGI", "Promoter-nonCGI")
#     annos <- c("genebody", "intergenic", "prom", "prom_cgi", "prom_noncgi")
#   }
# 
#   colors <- gg_color_hue(length(annos))
#   
#   data.frame(row.names = annos,
#              context = contexts,
#              color = gg_color_hue(length(annos)))
# }
# 
# 
# ############################ gg box plot by annotation - uses a df_anno data.frame
# gg_box <- function(dt,aesX,aesY, global=NULL,...){
#   
#   df_anno <- anno_df()
#   lvls <- rownames(df_anno)
#   lvls %<>% .[. %in% as.character(unique(dt$anno))]
#   
#   ## keep their order
#   df_anno %<>% .[lvls,]
#   
#   dt$anno %<>% factor(., levels = lvls, ordered = TRUE)
#   dt %<>% .[order(dt$anno, decreasing = TRUE),]
#   
#   p <- ggplot(dt) 
#   if(!is.null(global)){
#     p <- p + geom_rect(aes(xmin=min(dt[,aesX]), xmax=max(dt[,aesX]), ymin=global[1], ymax=global[6]), fill="grey80")
#     p <- p + geom_rect(aes(xmin=min(dt[,aesX]), xmax=max(dt[,aesX]), ymin=global[2], ymax=global[5]), fill="grey60")
#   }
#   
#   p <- p + theme_bw() +
#     geom_boxplot(aes_string(x=aesX, y=aesY,  col=aesX),alpha=0.7, show.legend = FALSE) +labs(...) +
#     geom_violin(aes_string(x=aesX, y=aesY, fill=aesX, col=aesX),alpha=0.7, show.legend = FALSE) +
#     scale_fill_manual(values = df_anno$color, label=df_anno$context, guide_legend(title="Context")) +
#     scale_color_manual(values = df_anno$color, label=df_anno$context, guide_legend(title="Context")) +
#     scale_x_discrete(labels=df_anno$context) + theme(axis.text.x = element_text(face="bold")) + coord_flip()
#   # scale_x_discrete(position = "top")
#   
#   return(p)
# }

############################ install someWrappers if not on my mac
if(!grepl('alabadi',Sys.info()['user'])){
  devtools::install_github("ajabadi/someWrappers")
}

############################ change chromosome format
chr_conv <- function(chr_vec, style = "UCSC"){ ## called by tri_stats; assign_strand_context; 
  ## convert chromosome names b/w UCSC and ensembl
  
  ## INPUT chr_vec: a vector of canonical + mt chromosomes only
  ## INPUT style: one of c("UCSC", "ensembl")
  ## OUTPUT converted char_vec
  
  if(style=="ensembl"){
    chr_vec <- ifelse(chr_vec=="chrM","MT",substring(chr_vec,first = 4))
  } else if(style=="UCSC"){
    chr_vec <- ifelse(chr_vec=="MT","chrM",paste0("chr",chr_vec))
  } else{stop("style must be one of c('UCSC', 'ensembl')")}
  
  return(chr_vec)
}


############################ assign strand and trinucleotide context to sample cov files
assign_strand_context <- function( dt = cov_file, ref_genome = BSgenome.Hsapiens.UCSC.hg38){
  ## assign strand and c context to calls using sequence info at loci
  ## INPUT dt: cov file which has 'start', 'end' (both equal) and 'chr' in ensemble format
  ##    chr    start      end 
  ##    10    10293020  10293020
  ## INPUT ref_genome: a BSgenome sequence
  ## OUTPUT updated dt withv$c_context column added
  ##    chr    start      end     c_context
  ##    10    10293020  10293020    CHH
  ## calculate the penta-nucleotide sequence at the call loci - 2 ahead for + and 2 before for -
  dt %<>% .[,tri_seq:=as.character(getSeq(x=ref_genome,
                                          names=chr_conv(chr, style = "UCSC"),
                                          start=start-2,
                                          end=end+2))]
  ## make sure the center one is a C or G
  if(!all(substring(dt$tri_seq,3,3) %in% c("C", "G"))) stop("some trinucleotide sequences at call loci do not match C or G")
  ## add strand based on whether ithe center one is C or not (G)
  dt[,strand:=ifelse(substring(tri_seq,3,3)=="C","+","-")]
  ## look down/up-stream for forward/reverse strand and reverse complement for reverse strand
  dt[,tri_seq:=ifelse(strand=="+", substring(tri_seq,3,5), 
                      as.character(Biostrings::reverseComplement(DNAStringSet(substring(tri_seq,1,3)))))]
  ## assign context
          ############################ get a trinuceotide character vector and assign context to it
          
          cytosine_context <- function(trinuc_vec){
            ##TODO
            
            if (!all(substring(trinuc_vec,0,1)=="C")) stop("not starting with C")
            out <- rep("CHH", length(trinuc_vec))
            out[startsWith(trinuc_vec,"CG")] <- "CpG"
            out[!startsWith(trinuc_vec,"CG") & endsWith(trinuc_vec,"G")] <- "CHG"
            
            return(out)
          }
          
  dt[,c_context:=cytosine_context(tri_seq)]
  dt[,c("strand", "tri_seq"):=NULL]
  return(dt)
}

############################ read bed files and add trinucleotide stats
read_bed <- function(bed_file=anno_file, genome = BSgenome.Hsapiens.UCSC.hg38, valid_chrs = opts$valid_chrs, subset=FALSE){
  ## read the annotation file and add cpg stats (number and density for a given id) to it
  
  ## INPUT bed_file: full path to a bed file of form (not necessarily with solumn names):
  ##    chr   start   end  strand  id   annotation
  ## INPUT valid_chrs: the chromosomes to keep, usually c(1:22, "MT") for hs
  ## INPUT genome: a BSgenome
  ## INPUT subset: FALSE, or an integer to subset the bed files up to, for test runs
  ## OUTPUT an annotation file with colnames: "id", "chr", "strand", "start", "end"
  
    out <- fread(bed_file ,sep="\t", header=F, verbose=F) %>%
    ## remove the  annotation
    .[,c(-6)] %>%
    ## rename the columns
    setnames(c("chr","start","end","strand","id"))
    
    if(subset){ ## must be an integer
      out %<>% .[1:subset]
    }
    ## add cpg_content  cpg_density  chg_content ... to data.table
    out <- tri_stats(anno_dt = out, genome = genome, chrs=valid_chrs)
  ## re-order
  # .[,c( "id", "chr", "strand", "cpg_total", "cpg_dens", "start", "end")]
  
  return(out)
}

############################ get annotation data.table and add cytsoine contexts content and density 
tri_stats <- function(anno_dt, genome=BSgenome.Hsapiens.UCSC.hg38, chrs = c(1:22, "MT")){
  ## 1. anno_dt!=NULL: gets a genomic annotation data.table and adds cytosines content/density for every CpG, CHG, CHH context
  ## 2. anno_dt=NULL: get a reference genome and output a data.table of cytosines content
  # cpg      chg      chh
  # 1: 12152703 17313475 67302871
  # 2:  9622606 13878778 55076340
  
  
  ## INPUT anno_dt: a data.table of at least these columns: (or NULL)
  #     chr strand  start     end
  # 1: chr1      + 923928  944581
  # 2: chr1      - 944204  959309
  # 3: chr1      + 960587  965715
  # 4: chr1      + 966497  975865
  # 5: chr1      - 975204  982093
  # 6: chr1      - 998962 1000172
  ## INPUT genome: reference genome
  ## INPUT: chrs: chromosomes to keep
  
  ## OUTPUT: 1. updated anno_dt with cytosine contexts (CpG, CHG, CHH)
  #     chr strand cpg_total   cpg_dens  start     end  cpg_content  cpg_density  chg_content ...
  # 1: chr1      +      1178 0.05703496 923928  944581    876           0.1 
  # 2: chr1      -      1237 0.03018668 944204  959309 
  # 3: chr1      +       371 0.07233379 960587  965715  
  # 4: chr1      +       360 0.03842459 966497  975865  
  # 5: chr1      -       602 0.05341074 975204  982093  
  # 6: chr1      -       176 0.12799339 998962 1000172  
  ## OUPUT 2. explained above
  
            ## used function
            ############################ get DNAStringSet and retunr a data.table of trincu content and density
            trinuc_content <- function(dnastring_set=seq){ ## called by tri_stats
              ## get a DNAStringSet and calculate cpg contents in each element
              
              ## INPUT dnastringset: a DNAStringSet object
              ## OUTPUT an integer of the number of cpgs
              
              trinucs <- trinucleotideFrequency(dnastring_set)
              trinucs %<>% (function(x) x[,startsWith(colnames(x),"C")])
              cpgs <- trinucs[,c("CGA","CGC","CGG","CGT")] %>% rowSums()
              chgs <- trinucs[,c("CAG","CCG","CTG")] %>% rowSums()
              chhs <- trinucs %>% .[ ,! colnames(.) %in% c(colnames(cpgs), colnames(chgs))] %>% rowSums()
              width_seq <- width(dnastring_set)
              return(data.table(cpg_content=cpgs, cpg_density=cpgs/width_seq, chg_content=chgs, chg_density=chgs/width_seq, 
                                chh_content=chhs, chh_density=chhs/width_seq))
            }

  
  if(!is.null(anno_dt)){
    
    ## keep only canonical and mt chromosomes
    anno_dt %<>%.[chr %in% chrs]
    seq <- with(anno_dt,BSgenome::getSeq(genome,names=chr_conv(chr, style = "UCSC"),start=as.integer(start),end=as.integer(end), strand=strand))
    
    anno_dt %<>% cbind(.,trinuc_content(seq))
    
    out <- anno_dt
  } else { ## if whole genome stats are required (anno_dt=NULL)
    ref_genome <- BSgenome::getSeq(genome,names = chr_conv(chrs, style = "UCSC"))
    
    genome_trinucs <- trinuc_content(ref_genome)
    genome_c_stats <- genome_trinucs[,c("cpg_content","chg_content", "chh_content")]
    out <- genome_c_stats
  }
  return(out)
}

############################ create arbitrary genomic windows

genomic_window <- function(genome=BSgenome.Hsapiens.UCSC.hg38,width=3000L,
                           step=c(600L,1500L,3000L), strand="+", valid_chr=c(1:22)){
  ## creates fixed length sliding/non-overlapping windows from a ref. genome
  ## if step==width,  it will be a non-overlapping one
  
  ## INPUT genome: Formal class 'BSgenome' [package "BSgenome"] 
  ## INPUT type: one/both of =c("sliding","tile")
  ## INPUT width: width of window
  ## INPUT step: for sliding window
  ## OUTPUT a list with non-overlappings in $sliding_step_STEP(s)
  
  out <- list()
  gr <- GRanges(seqinfo(BSgenome.Hsapiens.UCSC.hg38))
  
  for (step_i in step){
    out[[paste0("window_",width,"_",step_i)]] <- slidingWindows(gr, width = width, step = step_i) %>% as.data.table() %>% 
      .[,3:5] %>% .[,strand:=strand] %>% setnames(.,"seqnames","chr") %>% .[,chr:=substring(chr,4)] %>%.[chr%in%valid_chr,] %>% setkey(chr,start,end) 
  }
  
  return(out)
}
############################ weight calculator
weight_se <- function(r, total ){
  se2 <- r*(1-r)/(total)
  w <- 1/se2
  return(w)
}
############################ match enhancer and genes based on bed files and a window size
enhancer2geneMatch <- function(enhancer=enhancers, ## chr start end strand id
                               genebody=genebodies, ## chr start end strand id
                               size=1e6){ 
  ## returns an overlap table that you can dcast on gene_id~enhancer_id value.var='dist'
  
  gene_cenetred <- genebody %>% .[,`:=`(start=start-size, end=start+size, tss=start)]
  gene_cenetred[start<0,start:=0] %>% setkey(chr,start,end)
  ## long ones first, keyed, fat one second
  ovlap <- foverlaps(enhancer, gene_cenetred) ## overlaps with many, we want the closest
  # head(ovlap) ## we want the enhancer start to be close to be closest to 
  ## for the same overlaps (chr, start, end, i.start, i.end) we want the id with the least abs(tss-i.start)
  ovlap[,dist:=abs(tss-i.start)]
  ## sort by decreasing dist
  ovlap %<>% .[base::order(dist)] %>% setnames(c("id", "i.id"),c("gene_id", "enhancer_id"))
  ## get a lookup table
  match_all <- ovlap[,c("gene_id", "enhancer_id", "dist")]
  match_nearest <- match_all%>% .[!base::duplicated(., by=c("enhancer_id"))]
  return(list("match_all"=match_all, "match_nearest"=match_nearest))
}

#################### function to spot the well names using underscores b/w which the names can be found.
## INPUT:
## char_vec: characte vector of unprocessed cell names
## underscores: the index of the underscores. e.g. for cell A101 with raw name 'abc_zx_A101_ghf_dy', it is c(2,3)
## as the name is b/w the 2nd and the 3rd underscores

name_cells_by_wells <- function(char_vec,
                                underscores = c(1,2)){ ## index of undersores after and before whose the well name exists
  pos_underscore <- gregexpr('_', char_vec)
  ## add the zeroth underscore at 0
  pos_underscore = lapply(pos_underscore, function(x) c(0,x))
  underscores <- underscores+1 ## add one to adjust for addition of zerp
  substr_fun <- function(x,y, pos) substr(x, y[pos[1]]+1,y[pos[2]]-1)
  names <-  mapply(substr_fun,
                   x=char_vec,
                   y=pos_underscore, MoreArgs = list(pos=underscores))
  return(unname(names))
}

## in-line test

c01 <- c("A5_asas", "A9_shdghdg")
c12 <- c("x10_A61_hgh","ghsd_A5_hjhsd")

if(!all(name_cells_by_wells(c01, c(0,1))==c("A5","A9"), 
        name_cells_by_wells(c12, c(1,2))==c("A61", "A5"))){
  stop("name_cells_by_wells test failed")
}
rm(c01,c12)
