## function to calculate weighted mean and variance and return as rowData(sce)$wtd__var and $wtd_mean

wtd_var = function(sce=sce, x_assay = "rates", w_assay = "weights", cov = 0){
  
  ## filter by coverage
  sce = sce[rowSums(is.na(assay(sce, x_assay)))/dim(sce)[2]>= cov,]
  x = assay(sce, x_assay)
  w = assay(sce, w_assay)
  w_var = w_mean <- numeric()
  
  run.time = system.time({
    for (i in 1:dim(x)[1]){
      w_var[i] =  wtd.var(x = x[i,], weights = w[i,], na.rm = TRUE )
      w_mean[i] = wtd.mean(x = x[i,], weights = w[i,], na.rm = TRUE)
    }
  })
  ## NA variance is 
  rowData(sce)$wtd_var = w_var
  rowData(sce)$wtd_mean = w_mean
  sce = sce[order(rowData(sce)$wtd_var, decreasing = TRUE, na.last = TRUE),]
  return(sce)
}