library(biomaRt)
require(topGO)
require(org.Hs.eg.db)



# mouseMart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
#
# # mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# res <- getBM(attributes=c("ensembl_gene_id", "external_gene_name") , mart=mouseMart)
#
# mart <- res$external_gene_name %>% set_names(res$ensembl_gene_id)
# ##
comp=1
dat <- "X"
trms <- c('regul', 'morph')
# spls.obj ## has symbol gene names
# ##
# ## ID "ensembl" "symbol" "genename"
goEnrich <- function(spls.obj, dat, comp, mapping="org.Hs.eg.db", algorithm="classic", trms=c('regul', 'morph', 'differe', 'develop'), ID="ensembl", minEnrich=1.5, topNodes=20, boxfill=gg_color_hue(2)[2],
                     colMatch="purple", colNoMatch="black"){
  maxKS=10^(-minEnrich)
  ## a named vector of markers loadings
  markers <- selectVar(spls.obj, comp = comp)[[dat]]$value 
  markers <- markers$value.var %>% set_names(rownames(markers))
  ## feasible genes
  feas_genes <- spls.obj[[dat]] %>% colnames()
  
  
  
  selection <- function(allScore){ return(allScore != 0)} # function that returns TRUE/FALSE
  allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=feas_genes, mapping=mapping, ID=ID)
  GOdata <- new("topGOdata",
                ontology="BP",
                allGenes=markers,
                annot=annFUN.GO2genes,
                GO2genes=allGO2genes,
                geneSel=selection,
                nodeSize=10)
  
  ## use rank info, irrelevant here
  
  results.ks <- runTest(GOdata, algorithm=algorithm, statistic="ks")
  goEnrichment <- GenTable(GOdata, KS=results.ks, orderBy="KS", topNodes=topNodes)
  goEnrichment <- goEnrichment[goEnrichment$KS<maxKS,]
  goEnrichment <- goEnrichment[,c("GO.ID","Term","KS")]
  # goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
  # goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
  # goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
  goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
  goEnrichment$KS <- as.numeric(goEnrichment$KS)
  
  ## plot
  
  require(ggplot2)
  
  trm_match <- vector(length = dim(goEnrichment)[1])
  for (trm in trms){
    trm_match <- trm_match| grepl(trm, goEnrichment$Term)
  }
  trm_match <- rev(trm_match)
  
  trm_col <- ifelse(trm_match, colMatch, colNoMatch)
  # highlight <- rev(highlight)
  ggplot(goEnrichment, aes(x=Term, y=-log10(KS))) +
    stat_summary(geom = "bar", fun.y = mean, position = "dodge", fill=trm_col) +
    xlab("") +
    ylab("GO Term Enrichment") +
    scale_y_continuous(breaks = seq(0, max(-log10(goEnrichment$KS)), by = 0.5)) +
    theme_bw() +
    theme(
      legend.position='none',
      legend.background=element_rect(),
      plot.margin = margin(0, 0, 0, 0, "cm"),
      plot.title=element_text(angle=0, size=24, face="bold", vjust=1),
      axis.text.x=element_text(angle=0, size=12, face="bold", hjust=1.50),
      axis.text.y=element_text(angle=0, size=12, face="bold", vjust=0.5, color = trm_col),
      axis.title=element_text(size=12, face="bold", colour = "grey20"),
      legend.key=element_blank(),     #removes the border
      legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
      legend.text=element_text(size=18),  #Text size
      title=element_text(size=18)) +
    guides(colour=guide_legend(override.aes=list(size=2.5))) +
    coord_flip()
  # return(p)
}

