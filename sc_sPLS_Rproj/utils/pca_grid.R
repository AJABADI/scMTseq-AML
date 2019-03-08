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
pca_grid = function(pca_obj=pca_res, col='well', diff_samples = NULL,
                    default_label = "Matching", diff_label = "RNA-seq only", point_label=F, top = ""){
  
  gg_pca = function(pca=pca_res, PCs = c(1,2), label = point_label ){
    ## variates
    mat = as.data.frame(pca$variates$X)
    ## add well row
    mat$well = factor(substr(rownames(mat),0,1))
    
    if(is.null(diff_samples)){
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
        if(col!='well'){
          p = p + geom_point(aes(x=mat[,PCs[1]], y=mat[,PCs[2]]), col =col)
        } else{
          p = p + geom_point(aes(x=mat[,PCs[1]], y=mat[,PCs[2]], col = well))+
            scale_colour_manual(values = myDistinctColors)
        }
        
      }
    } else {
      mat$colors = default_label
      mat[diff_samples,]$colors = diff_label
      p = ggplot(mat) + geom_point(aes(x=mat[,PCs[1]], y=mat[,PCs[2]], col = colors))  +
        labs(x = paste0('PC',PCs[1]), y= paste0('PC',PCs[2]))+
        theme_bw() + 
        theme(plot.title = element_text(hjust = 0.5)) + 
        guides(col=guide_legend(title="cells"))
    }
    return(p)
  }
  
  
  exp_var = data.frame(pca_obj$explained_variance)
  p1 = ggplot(exp_var) +geom_col(aes(x=rownames(exp_var), y=exp_var[,1]))+
    labs(x='PCs', y= 'Explained Variance') + ggtitle('Explained Variance') +
    theme(plot.title = element_text(hjust = 0.5))
  
  extract_legend <-function(a.gplot=p2){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  ## legend for all plots
  p2 = gg_pca(pca=pca_obj, PCs = c(1,2))
  all_legend = extract_legend(p2)
  
  p2 = p2+ theme(legend.position="none")
  p3 = gg_pca(pca=pca_obj, PCs = c(1,3))+ theme(legend.position="none")
  p4 = gg_pca(pca=pca_obj, PCs = c(2,3))+ theme(legend.position="none")
  
  grid = grid.arrange(arrangeGrob(p1,p2,p3,p4, ncol = 2), all_legend,nrow=1,widths = c(11,1), top =textGrob(top,gp=gpar(fontsize=16,font=3, col = "blue")))
  return(grid)
}
