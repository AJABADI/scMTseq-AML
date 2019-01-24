## function to plot sample PCA plots and color them by well group
## so that possible tehnical effects can show

pca.well = function(pca.obj = pca.res, PCs= c(1,2)){
  
  ## PCs from pca object
  df = as.data.frame(pca.res$variates$X)
  ## row to be the first letter of sample name (e.g. "A" for "A12" and )
  df$row = factor(substr(rownames(df),0,1))
  ################################# define custom colors
  ## 14 distinct colors
  
  myDistinctColors = c("#000000", ## black
                       "#3cb44b", ## green
                       "#000075", ## marine blue
                       "#42d4f4", ## cyan
                       "#f58231", ## orange
                       "#e6194B", ## red but pink
                       "#ffe119", ## yellow
                       "#911eb4", ## purple
                       "#9A6324", ## brown
                       "#808000", ## olive
                       "#f032e6", ## magneta (purple-pink)
                       "#800000", ## maroon
                       "#a9a9a9", ## grey
                       "#bfef45"  ## lime
                       
  )
  #################################
  p = ggplot(df, aes(df[,PCs[1]], df[,PCs[2]], col=row)) +
    geom_point(size=3)+
    geom_text(aes(label=rownames(df), col=row), hjust=0.2, vjust=-0.6) +
    scale_colour_manual(values = myDistinctColors) +
    labs(x = paste0("PC ", PCs[1]," : ",
                    round(pca.res$explained_variance[PCs[1]]*100, digits = 1),"%"),
         y = paste0("PC ", PCs[2]," : ",
                    round(pca.res$explained_variance[PCs[2]]*100, digits = 1),"%"),
         title = "PCA of scRNAseq cells coloured by well rows") +
    guides(col=guide_legend("well row"))
  return(p)
}