## function to create sce from all annotations in sce and
## also compute the coverage of cells by each gene in rowData
## assumes 'id' is in the first column

met_sce = function(met=met, annot="genebody", tot_c = NULL, cov = 0.1){
  rates = dcast(met[anno==annot], id~sample, value.var = 'rate')
  weights = dcast(met[anno==annot], id~sample, value.var = 'weight')
  ## modify this when you analyse windows - total number of possible methylation calls
  if(!is.null(tot_c)){
    c_content = dcast(met[anno==annot], id~sample, value.var = tot_c)
  }
  ## add id as rowname
  rates %<>% as.data.frame(., row.names = .[,1]) %>% .[,-1]
  weights %<>% as.data.frame(., row.names = .[,1]) %>% .[,-1]
  ## create sce
  sce = SingleCellExperiment(assays = SimpleList(rates = rates,
                                                 weights = weights))
  ## add cell coverage as rowData
  rowData(sce)$cov = rowSums(!is.na(rates))/ncol(rates)
  
  ## add weighted variance and mean as rowData
  sce %<>% wtd_var(., x_assay = "rates", w_assay = "weights", cov = cov)
  
  return(sce)
}
  
  
