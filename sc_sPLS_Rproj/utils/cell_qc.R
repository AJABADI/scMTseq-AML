## plot cells library size vs number of genes expressed and filter cells by percentage of genes expressing
## 
cell_qc = function(sce = rna_sce, pct_genes_in_cells = 10){
  stopifnot(pct_genes_in_cells>0 & pct_genes_in_cells<=100)
  df_cell_qc = data.frame(row.names = colnames(sce),
    tpm = sce$total_counts/1e6,
    genes_expressed = sce$total_features_by_counts/1e3
  )
  
  ## plot of library size vs number of detected genes
  p=ggplot(df_cell_qc, aes(x = tpm, y = genes_expressed)) +
    geom_point(size = 3, col=color.mixo(1)) + 
    labs(x = 'Library Size (Million)' , y=' Number of genes expressed (x 1000)', title = '# of genes expressed against library size for all cells.') +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_label_repel(
      aes(label = rownames(df_cell_qc)),
      col = 'darkorange',
      fontface = 'bold',
      box.padding = unit(0.35, 'lines'),
      point.padding = unit(0.5, 'lines'),
      segment.color = 'grey50'
    )
  ## histogram of percentage of detected genes for all cells
  hist_pct_genes = sce$pct_total_features_by_counts %>% as.data.frame() %>% 
    ggplot() + geom_histogram(aes(.), bins=15, fill='grey70', col='blue') +
    labs(x='Percent of detected genes')
  
  ## if the minimum detection limit is less than the threshold, filter
  if(pct_genes_in_cells>min(sce$pct_total_features_by_counts)){
    y_intercept = (dim(sce)[1]*pct_genes_in_cells/100)/1e3
    p = p+ geom_hline(yintercept = y_intercept, col='red')
    cells_pass = colnames(sce[,sce$pct_total_features_by_counts>=pct_genes_in_cells])
    sce = sce[,cells_pass]
  }
  
  return(list('sce'=sce,'plot' = p, 'hist'=hist_pct_genes))
}