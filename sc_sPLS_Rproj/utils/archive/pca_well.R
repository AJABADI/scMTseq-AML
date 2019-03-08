## function to plot sample PCA plots and color them by well group
## so that possible tehnical effects can show

pca_well = function(pca_obj = pca_res, PCs= c(1,2)){
  
  ## PCs from pca object
  df = as.data.frame(pca_res$variates$X)
  ## row to be the first letter of sample name (e.g. "A" for "A12" and )
  df$row = factor(substr(rownames(df),0,1))

  #################################
  p = ggplot(df, aes(df[,PCs[1]], df[,PCs[2]], col=row)) +
    geom_point(size=3)+
    geom_text(aes(label=rownames(df), col=row), hjust=0.2, vjust=-0.6) +
    scale_colour_manual(values = myDistinctColors) +
    labs(x = paste0("PC ", PCs[1]," : ",
                    round(pca_res$explained_variance[PCs[1]]*100, digits = 1),"%"),
         y = paste0("PC ", PCs[2]," : ",
                    round(pca_res$explained_variance[PCs[2]]*100, digits = 1),"%"),
         title = "PCA of scRNAseq cells coloured by well rows") +
    guides(col=guide_legend("well row"))
  return(p)
}