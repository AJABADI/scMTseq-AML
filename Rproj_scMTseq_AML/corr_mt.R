## function to calculate correlation of rna and methylation
## correlation of each gene and p-value is stored in rowData(sce)$corr_rna and $pvalue
## correlation are calculated using top top_prop proportion of common genes
## based on v_hat
corr_mt = function(rna = rna_match, met = sce_prom, top_prop =  0.3){
  
  expr = logcounts(rna)
  rates = assay(met, "rates")
  Weights = assay(met, "weights")
  
  ## intersect genes
  intersect = Reduce(intersect , list(rownames(expr), rownames(rates)))
  ## if top hvgs required - filter using top proportion of common genes
  if(!is.null(top_prop)){intersect %<>% .[1:floor(length(.)*top_prop)]}
  expr = expr[intersect,]
  rates = rates[intersect,]
  Weights = Weights[intersect,]
  
  pearson_genewise = data.frame(matrix(ncol = 4, nrow = dim(rates)[1]),
                                row.names = row.names(rates))
  for (i in intersect){
    X = expr[i,];  Y = rates[i,]; W = Weights[i,]
    ## remove NA cells
    non_NA_cells = names(W[which(!is.na(W))])
    X = X[non_NA_cells];  Y = Y[non_NA_cells]; W = W[non_NA_cells]
    X %<>% as.numeric()
    Y %<>% as.numeric()
    W %<>% as.numeric()
    ## if var(Y) or Var(X) is zero there is no correlation
    pearson_genewise[i,] = ifelse(var(Y)&var(X),weights::wtd.cor(X,Y,W),0)
  }
  ## initialise rowData
  rowData(met)$corr_rna = NA
  rowData(met)$corr_pvalue = NA
  ## replace computed ones
  rowData(met[intersect,])$corr_rna = pearson_genewise[,1]
  rowData(met[intersect,])$corr_pvalue = pearson_genewise[,4]
  
  ## corr of variances
  ## using total variance
  corr_var_tot = weights::wtd.cor(rowData(rna[intersect,])$total,
                                  rowData(met[intersect,])$corr_rna,
                                  rowData(met[intersect,])$cov)
  ## using biological variance
  corr_var_bio = weights::wtd.cor(rowData(rna[intersect,])$bio,
                                  rowData(met[intersect,])$corr_rna,
                                  rowData(met[intersect,])$cov)
  ## add to metadata
  metadata(met)$corr_var_tot_wtd_by_cov = corr_var_tot
  metadata(met)$corr_var_bio_wtd_by_cov = corr_var_bio
  return(met)
}