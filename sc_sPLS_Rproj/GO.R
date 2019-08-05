## enrichment
# library('biomaRt')
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# require(topGO)
# require(org.Hs.eg.db)

for (anno in names(list_geneset$canonical$fullset)[1:4]){
  for (view in c('transcriptome', 'methylome')){
    anno='genebody'
    view='transcriptome'
## supply gene names
genes <- list_geneset$canonical$fullset[[anno]][[view]]
genes <- genes[!duplicated(genes)]
all_genes <- rownames(transcriptome)
if(view!='transcriptome') all_genes <- unique(met_exp[anno==anno]$id)
  all_genes <- genes
genes <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","description"),values=genes,mart= mart)
all_genes <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","description"),values=all_genes,mart= mart)
saveRDS(all_genes, file=paste0('all_genes-',anno,view,'.Rds'))
## a p=value
genes_id <- rep(0.049, length(genes$ensembl_gene_id))
names(genes_id) <- genes$hgnc_symbol


selection <- function(allScore){ return(allScore < 0.05)} # function that returns TRUE/FALSE for p-values<0.05
allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=all_genes$ensembl_gene_id, mapping="org.Hs.eg.db", ID="symbol")
GOdata <- new("topGOdata",
              ontology="BP",
              allGenes=genes_id,
              annot=annFUN.GO2genes,
              GO2genes=allGO2genes,
              geneSel=selection,
              nodeSize=10)

## use rank info, irrelevant here

results.ks <- runTest(GOdata, algorithm="classic", statistic="ks")
goEnrichment <- GenTable(GOdata, KS=results.ks, orderBy="KS", topNodes=20)
goEnrichment <- goEnrichment[goEnrichment$KS<0.05,]
goEnrichment <- goEnrichment[,c("GO.ID","Term","KS")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
goEnrichment$KS <- as.numeric(goEnrichment$KS)

## plot
pdf(paste0('enrichment-',anno,view,'.pdf'))
require(ggplot2)
ggplot(goEnrichment, aes(x=Term, y=-log10(KS))) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("Enrichment") +
  ggtitle(paste0(anno,'-' ,view)) +
  scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$KS)), by = 2), 1)) +
  theme_bw(base_size=24) +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=24, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=18, face="bold", hjust=1.10),
    axis.text.y=element_text(angle=0, size=18, face="bold", vjust=0.5),
    axis.title=element_text(size=24, face="bold"),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=18),  #Text size
    title=element_text(size=18)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()
dev.off()
  }
}