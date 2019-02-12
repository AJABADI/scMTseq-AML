## function to calculate weighted mean and variance and return as rowData(sce)$wtd__var and $wtd_mean

wtd_var = function(sce=sce, x_assay = "rates", w_assay = "weights", cov=NULL){
  
  if(!is.null(cov)){ ## if filtering by coverage desired
    ## add cell coverage as rowData
    rowData(sce)$cov = rowSums(!is.na(rates))/ncol(rates)
    
    ## filter by coverage
    sce %<>% .[rowData(.)$cov>=cov,]
  }
  x = assay(sce, x_assay)
  w = assay(sce, w_assay)
  w_var = w_mean <- numeric()
  
  run.time = system.time({
    for (i in 1:dim(x)[1]){
      w_var[i] =  wtd.var(x = x[i,], weights = w[i,], na.rm = TRUE )
      w_mean[i] = wtd.mean(x = x[i,], weights = w[i,], na.rm = TRUE)
    }
  })
  ## add to rowData
  rowData(sce)$wtd_var = w_var
  rowData(sce)$wtd_mean = w_mean
  sce = sce[order(rowData(sce)$wtd_var, decreasing = TRUE, na.last = TRUE),]
  return(sce)
}