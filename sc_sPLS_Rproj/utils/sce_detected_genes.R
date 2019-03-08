## function to filter undected genes and calculate the percentage of total features detected in each cell
## and add to colData

sce_detected_genes = function(sce) {
  sce = sce[rowData(sce)$mean_counts>0,]
  sce$pct_total_features_by_counts = 100*sce$total_features_by_counts/dim(sce)[1]
  return(sce)
}