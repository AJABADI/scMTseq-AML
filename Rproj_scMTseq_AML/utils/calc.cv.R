#################### function to create expression summaries (mean, sd, cv)
## INPUT:
      ## count.data: a matrix of counts
      ## mean.expr.cutoff: below which HVGs can't be defined

library(ggplot2); theme_set(theme_bw())
calc.cv <- function(count.data=logcounts(sce.rna.scran),
                    mean.expr.cutoff = 3, ## put cv = NA for less
                    mean.quant = 0.25, ## mean expression qunatile to deflate CV
                    n.hvgs = NULL, ## number of top-cv genes to label HVG
                    type = 'Expression' ## c('Expression, 'Methylation')
){
  stopifnot(is.numeric(n.hvgs)||!is.null(n.hvgs))
  # count.data[is.na(count.data)]=0
  ## create a dataframe of mean expression and SD
  df <- data.frame(mean.expr=rowMeans(count.data, na.rm = T), SD = apply(count.data,1,sd, na.rm=T))
  ## calcuate CV
  df$CV = df$SD/(df$mean.expr + quantile(df$mean.expr,mean.quant ))
  df$CV[df$mean.expr<mean.expr.cutoff] = NA
  df = df[order(df$CV, decreasing = T),]
  if(!is.null(n.hvgs)){
    df$HVG = 'NO'
    df$HVG[df$mean.expr>mean.expr.cutoff][1:n.hvgs] = 'YES'
  }
  
  ## mean expression vs CV
  xtitle = paste0('Mean ', type)
  Title = paste0(xtitle, ' vs CV')
    if(!is.null(n.hvgs)){
      p<- ggplot(df) + geom_point(aes(x=mean.expr, y=CV, col=HVG))+
        labs(x = xtitle , y='Coefficient of Variation', title=Title) +
        guides(colour=guide_legend(title = 'HVG'))
    } else {
      p<- ggplot(df) + geom_point(aes(x=mean.expr, y=CV), col = 'blue')+
        labs(x = xtitle , y='Coefficient of Variation', title=Title)
    }
  
  return(list('expr.smry' = df, 'plot.cv' = p))
}