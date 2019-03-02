#################### functions to plot pca objects using ggplot2
## INPUTS: 
## pca.obj: a mixOmics' pca object with at least 3 components
## top: a title
## OUTPUT:
## a grid plot of explained variance in each pcs and pcs against each other
library(gridExtra)
library(ggplot2)

# pca.obj=pca.res
# col='darkblue'
# diff.samples = NULL
# default.label = "Matching"
# diff.label = "RNA-seq only"
# point.label=F


# pca = pca.obj
# PCs = c(1,2)

## multi pca + variance plot
pca.grid = function(pca.obj=pca.res, col='darkblue', diff.samples = NULL,
                    default.label = "Matching", diff.label = "RNA-seq only", point.label=F, top = ""){
  
  gg.pca = function(pca=pca.res, PCs = c(1,2), label = point.label ){
    mat = as.data.frame(pca$variates$X)
    if(is.null(diff.samples)){
      p = ggplot(mat) +
        labs(x = paste0('PC',PCs[1]), y= paste0('PC',PCs[2]))+
        theme_bw() + 
        theme(plot.title = element_text(hjust = 0.5))
      
      if(label){
        well = substr(rownames(mat),0,1)
        p = p + 
          geom_text(aes(x=mat[,PCs[1]], y=mat[,PCs[2]],label=rownames(mat),
                        col = well),hjust=0, vjust=0) +
          scale_color_manual(values = color.mixo(1:length(well)))
      } else {
        p = p + geom_point(aes(x=mat[,PCs[1]], y=mat[,PCs[2]]), col =col)
      }
    } else {
      mat$colors = default.label
      mat[diff.samples,]$colors = diff.label
      p = ggplot(mat) + geom_point(aes(x=mat[,PCs[1]], y=mat[,PCs[2]], col = colors))  +
        labs(x = paste0('PC',PCs[1]), y= paste0('PC',PCs[2]))+
        theme_bw() + 
        theme(plot.title = element_text(hjust = 0.5)) + 
        guides(col=guide_legend(title="cells"))
    }
    return(p)
  }
  
  
  exp.var = data.frame(pca.obj$explained_variance)
  p1 = ggplot(exp.var) +geom_col(aes(x=rownames(exp.var), y=exp.var[,1]))+
    labs(x='PCs', y= 'Explained Variance') + ggtitle('Explained Variance') +
    theme(plot.title = element_text(hjust = 0.5))
  
  p2 = gg.pca(pca=pca.obj, PCs = c(1,2))
  p3 = gg.pca(pca=pca.obj, PCs = c(2,3))
  p4 = gg.pca(pca=pca.obj, PCs = c(1,3))
  grid = grid.arrange(p1,p2,p3,p4, top =textGrob(top,gp=gpar(fontsize=16,font=3, col = "blue")))
  return(grid)
}
