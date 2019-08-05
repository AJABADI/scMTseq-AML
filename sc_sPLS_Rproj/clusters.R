######## scRNAseq
library(data.table)
library(SingleCellExperiment)
sce_rna <- readRDS('/Users/alabadi/Projects/analysis/R/multiOmics/sc_sPLS_aml/output/parsed/rna_sce_normalised_red_dims.Rds')
met_exp <- readRDS('/Users/alabadi/Projects/analysis/R/multiOmics/sc_sPLS_aml/sc_sPLS_Rproj/02_Methylation/02_QC_Methylation/met_exp.Rds')
id_intersect <- Reduce(intersect, list(rownames(sce_rna), unique(met_exp$id)))
length(id_intersect)
####### assess the pairwise correlations for common highly variable genes
## sort sce_rna by hvgs
sce_rna %<>% .[rowData(.) %>% data.frame() %>% .[order(-.$bio),] %>% rownames(),]

## keep common ones
sce_rna %<>% .[rownames(.) %in% id_intersect] ## from 9394 to 8214
met_exp %<>% .[id %in% id_intersect] ## from 589123 to 62496

## calculate the pairwise correlation for common hvgs
hvg_corr <- function(sce_rna, met_exp, init_hvg=1000, anno='genebody',...){
  met_anno <- met_exp[anno==anno]
  met_anno %<>% .[order(-lci)]
  unique_met_anno <- met_anno[!duplicated(id)]
  met_hvgs <- unique_met_anno$id[1:init_hvg]
  rna_hvgs <- rownames(sce_rna)[1:init_hvg]
  common_hvgs <- intersect(met_hvgs, rna_hvgs)
  message(length(common_hvgs))
  return(with(met_anno[id %in% common_hvgs], stats::cor(bio, lci, use='pairwise.complete',...)))
}

met_exp %<>% .[order(-lci)]
unique_metexp <- met_exp[!duplicated(id)]
library(propr)
library(ggplot2)

library(gplots)
## create a propr class from logcounts
pr <- propr(logcounts(sce_rna))
M <- getMatrix(pr)
h <- heatmap(M, keep.dendro = TRUE)
row.clusters = as.hclust( h$Rowv )
n_clust <- 4
clusters <-factor(cutree(row.clusters,k=n_clust), levels = 1:n_clust)
colz <- clusters
levels(colz) <- gg_color_hue(n_clust)

heatmap(M,  dendrogram = 'row', trace='none')
## pca
ggplot(as.data.frame(reducedDim(sce_rna, 'pca_all'))) + 
  geom_point(aes(x=PC1, y=PC2, col=factor(clusters))) + guides(col=guide_legend(title = 'Cluster'))

## umap
# library(umap)
# k <- 4
# umap.res <- umap(d = t(logcounts(sce_rna)), method = 'naive', n_components=k)
# colnames(umap.res$layout) <- paste0('UMAP_',1:k)
# df <- as.data.frame(umap.res$layout)
# ggplot(df) + geom_point(aes(x=UMAP_1, y=UMAP_2 , col=factor(clusters[rownames(df)])))

## function to case methylation data from the data.table

dcast_met <- function(met_dt=met_exp, annotation='genebody', 
                      min.cov=10, ## min cell coverage
                      filt_sample=NULL, ## samples to filter
                      valueVar='rhat', ## value to cast
                      hvgs=FALSE){ ## FALSE or integer
  met_dt %<>% .[anno==annotation & cell_cov>=min.cov & !sample %in% filt_sample]
  
  if(hvgs){
    met_dt %<>% .[order(-lci)]
    unique_dt <- met_dt[!duplicated(id)]
    if(dim(unique_dt)[1]<hvgs){
      message('Number of hvgs enforced by dimensions')
      hvgs <- dim(unique_dt)[1]
    }
    hv_ids <- unique_dt$id[1:hvgs]
    met_dt %<>% .[id %in% hv_ids]
  }
  
  dcast(met_dt,id~sample, 
        value.var = valueVar ) %>% setkey(id) %>% as.matrix(rownames=TRUE)
}


met_exp <- readRDS("/Users/alabadi/Documents/_Projects/sc_sPLS/sc_sPLS_Rproj/02_Methylation/02_QC_Methylation/met_exp.Rds")
######## BSseq
## weighted euclidean norm distance matrix
wt_euc_dist <- function(met_exp, anno, minCov, hvgs=FALSE) {
  mat <- dcast_met(met_dt = met_exp, annotation = anno, min.cov = minCov,hvgs=hvgs, valueVar = 'rhat')
  wt <- dcast_met(met_dt = met_exp, annotation = anno, min.cov = minCov, hvgs=hvgs, valueVar = 'wij')
  calls <- dcast_met(met_dt = met_exp, annotation = anno, min.cov = minCov, hvgs=hvgs,  valueVar = 'calls')

  mat <- mat[,order(colnames(mat))]
  iter <- expand.grid(list(rownames(mat),colnames(mat), colnames(mat)))
  names(iter) <- c('site','cell1', 'cell2')
  iter %<>% .[with(iter,cell1!=cell2),]
  iter %<>% .[order(iter$cell1, rev(iter$cell2)),] 
  iter %<>% .[1:floor(dim(iter)[1]/2),]
  iter$wijj <- mapply(function(site,cell1,cell2) (sqrt(prod(wt[site,cell1],wt[site,cell2], na.rm = TRUE))), 
                      iter[,'site'], iter[,'cell1'], iter[,'cell2'])
  iter$r1 <-  mapply(function(site,cell1) (mat[site,cell1]), 
                      iter[,'site'], iter[,'cell1'])
  iter$r2 <-  mapply(function(site,cell2) (mat[site,cell2]), 
                      iter[,'site'], iter[,'cell2'])
  iter$calls <- mapply(function(site) (sum(calls[site,], na.rm = TRUE)),  iter[,1])

  iter <- as.data.table(iter, keep.rownames = F)
  iter[,sumwt:=sum(wijj, na.rm = TRUE), by=site]
  iter[,wijj:=wijj*calls/sum(wijj), by=site]
  iter[,l2_wt:=wijj*(r1-r2)^2]
  dists <- iter[, .(wt_euc_d=sqrt(sum(l2_wt, na.rm = TRUE))), by=c('cell1', 'cell2')]
  dists %<>% .[order(dists$cell1, dists$cell2)]
  dist_mat <- matrix(0, nrow = dim(mat)[2], ncol = dim(mat)[2])
  colnames(dist_mat) <- rownames(dist_mat) <- colnames(mat)
  # dist_mat[upper.tri(dist_mat)] <- dists$wt_euc_d
  dist_mat[lower.tri(dist_mat)] <- dists$wt_euc_d
  dist_mat <- dist_mat + t(dist_mat)
  return(dist_mat)
}

## MDS
## test

anno_cov <- rep(20,5) #c(20,20,20,20,20)
anno_hvgs <- rep(300,5) #c(300,300,300,300,300)
anno_context <- c('Genebody', 'Promoter - CGI', 'Promoter', "Promoter non-CGI", 'Intergenic')
anno_all <- c("genebody", "prom_cgi", "prom", "prom_noncgi", "intergenic")
anno_df <- data.frame(cov=anno_cov, context=anno_context, hvgs=anno_hvgs, row.names = anno_all )

dmat_list_hvg <- list()

for(i in rownames(anno_df)){
  message(i)
  dmat_list_hvg[[i]] <- wt_euc_dist(met_exp = met_exp, anno=i, minCov = anno_df[i,]$cov, hvgs = anno_df[i,]$hvgs)
}
saveRDS(dmat_list_hvg, file = 'dmat_list_hvg_300_cov_20.Rds')
file_backup('dmat_list_hvg_300_cov_20.Rds', where = NA, run_spec = 'distance matrices with 300 highly variable genes based on lci')

## get a list of cluster memberships for each annotation given a distance matrix

## plot the heatmaps and decide on clusters
for (i in rownames(anno_df)){
  titl <- paste0(anno_df[i,]$context, ':\n Hierarchical clustering of context methylome \ncoloured by transcriptome clusters')
  gplots::heatmap.2(dmat_list_hvg[[i]],  RowSideColors = as.character(colz),
                    margins = c(1,1), dendrogram ='row',trace='none' ,main = titl)
}

met_clust <- function(anno_df, dmat=dmat_list_hvg, n_clust){
  

met_clusters <- list()
for(i in rownames(anno_df)){
  h <- heatmap(dmat[[i]], keep.dendro = TRUE)
  row.clusters = as.hclust( h$Rowv )
  met_clusters[[i]] <- cutree(row.clusters,k=n_clust)
}


gplots::heatmap.2(dmat_list$intergenic,  RowSideColors = as.character(colz),margins = c(1,1), dendrogram ='row',trace='none' ,main = 'Hierarchical clustering of genebody methylome \ncoloured by transcriptome clusters')


mds_list <- list()
for(i in names(anno_cov)){
  mds <- cmdscale(dmat_list[[i]], k=5)
  mds_list[[i]]<- data.frame(mds) %>% set_colnames(paste0("MDS-",1:dim(mds)[2]))
  rm(mds)
}


mds_plot <- function(mds, anno, clusters, anno_df, comps=c(1,2)){
  xy = paste0('`MDS-',comps,'`')
  p <- ggplot(mds, aes_string(x=xy[1], y=xy[2])) + geom_point(aes(col=factor(clusters[rownames(mds)])), size=3) + 
    ggtitle(anno_df[anno,]$context) +  geom_text(aes(label=rownames(mds)),col='grey50', vjust=-0.7, hjust=0.7, size=3) + theme_bw() +
    guides(col=guide_legend(title = 'Transcriptome \nCluster')) + theme(plot.title = element_text(hjust=0.5))
  return(p)
}
}






j=1
for(i in names(anno_cov)){
  i <- names(anno_cov)[j]
  j=j+1
  mds <- mds_list[[i]]
  pdf(paste0('./03_Integration/img/mds_',i,'.pdf'))
  mds_plot(mds=mds,anno = i, clusters = clusters,anno_df = anno_df, comps = c(1,2))
  dev.off()
  # rm(mds)
}

anno='intergenic'
comps=c(1,2)
mds_plot(mds=mds_list[[anno]],anno = anno, clusters = clusters,anno_df = anno_df, comps = comps)


# ## seurat clustering
# library(Seurat)
# seurat <- CreateSeuratObject(raw.data = counts(sce_rna), project = "aml", min.cells = 0, min.features = 0)
# seurat <- NormalizeData(object = seurat, normalization.method = "LogNormalize", scale.factor = 10000)
# seurat <- FindVariableGenes(object = seurat, top.genes = 500)
# top200 <- head(x = VariableGene(object = seurat), 200)
# VariableGenePlot(seurat)
# seurat@var.genes
# seurat %<>% ScaleData() 
# seurat %<>% RunPCA()
# seurat <- FindClusters(object = seurat, resolution = 1.2)
# seurat@ident
# ggplot(as.data.frame(reducedDim(sce_rna, 'pca_all'))) +geom_point(aes(x=PC1, y=PC2, col=factor(seurat@ident[colnames(sce_rna)])))
# 
# ggplot(df) + geom_point(aes(x=UMAP_1, y=UMAP_2 , col=factor(seurat@ident[rownames(df)])))
# seurat %<>% RunUMAP(dims=1:5)
# markers <- FindMarkers(seurat, ident.1 = 0)
# markers

