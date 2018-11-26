#################### functions to plot pca objects using ggplot2
## INPUTS: 
      ## pca.obj: a mixOmics' pca object with at least 3 components
      ## top: a title
## OUTPUT:
      ## a grid plot of explained variance in each pcs and pcs against each other

gg.pca = function(pca=pca.scran, PCs = c(1,2), col='darkblue'){
  mat = as.data.frame(pca$variates$X)
  p = ggplot(mat) + geom_point(aes(x=mat[,PCs[1]], y=mat[,PCs[2]]), col =col)  +
    labs(x = paste0('PC',PCs[1]), y= paste0('PC',PCs[2]))+
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

## multi pca + variance plot
plot.pca = function(gg.pca=gg.pca, pca.obj=pca.scran, top = 'scran-normalised'){
  exp.var = data.frame(pca.obj$explained_variance)
  p1 = ggplot(exp.var) +geom_col(aes(x=rownames(exp.var), y=exp.var[,1]))+
    labs(x='PCs', y= 'Explained Variance') + ggtitle('Explained Variance') +
    theme(plot.title = element_text(hjust = 0.5))
  
  p2 = gg.pca(pca=pca.obj, PCs = c(1,2))
  p3 = gg.pca(pca=pca.obj, PCs = c(2,3))
  p4 = gg.pca(pca=pca.obj, PCs = c(1,3))
  grid = grid.arrange(p1,p2,p3,p4, top= textGrob(top, gp=gpar(col="darkorange", fontsize=18 )))
  return(grid)
}