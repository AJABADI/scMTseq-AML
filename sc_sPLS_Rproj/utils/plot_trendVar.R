## plot variance vs mean expression
plot_trendvar = function(sce, caption='', low='grey60', high='blue', alpha=0.75){
  row_data =  sce %>% rowData() %>% as.data.frame() ## a data.frame of rowData
  ## total variance (fold change compared to mean) vs log mean
  p = ggplot(row_data) + geom_point(aes(x= mean, y=total, col = bio), show.legend = TRUE, alpha=alpha) +
    geom_line(aes(x=mean, y = trendVar, fill = 'Fitted trend'),  size=2) +
    scale_fill_manual('', values = c(6), guide=guide_legend(override.aes = list(color=c('black')))) +
    scale_color_gradient(low=low, high = high, name = 'Bio. Variance')+
    labs(x='mean log expression', y = 'total variance', title = paste0('total variance vs log mean expression', caption))
  return(p)
}