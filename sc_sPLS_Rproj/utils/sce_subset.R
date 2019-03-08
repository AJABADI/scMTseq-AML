## function to subset sce to n_top features based on a row_data
sce_subset <- function(sce = rna_sce, row_data = 'total', n_top =1000){
  sce <- sce[unlist(order(rowData(sce)[row_data], decreasing = TRUE)),]
  sce <- sce[1:n_top,]
  return(sce)
}