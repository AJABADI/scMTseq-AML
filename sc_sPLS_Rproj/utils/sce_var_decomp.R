## function to decompose variance for sce and return results as rowData(sce)
## ordering genes by biological variance
sce_var_decomp <-  function(sce=rna_sce){
  fit <- trendVar(sce, method='loess', use.spikes=FALSE)  ## fit a mean-dependent loess to variance
  decomp <- as.data.frame(decomposeVar(sce, fit)) ## decompose variance to technical and biological
  decomp$trendVar <- fit$trend(fit$means) ## add the fitted variance column
  rowData(sce) %<>%  cbind(.,decomp) ## add to rowData(sce)
  
  return(sce)
}