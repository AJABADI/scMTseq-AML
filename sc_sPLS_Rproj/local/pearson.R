library(data.table)
library(stats)
tp <- system.time({
  pearson_genewise <- list()
  for (annot in annotations){
    message(annot)
    dt <- met_exp[id %in% common_genes,][anno==annot]
    
    dt %<>% .[cpg_cov>cov_cutoff]
    expr<- dcast.data.table(dt,id~sample, value.var = 'expression' )
    methyl <- dcast.data.table(dt,id~sample, value.var = 'rhat' )
    ## TODO add weights here
    suppressWarnings( pearson_genewise[[annot]] <- sapply(common_genes, FUN = function(x) cor(as.numeric(expr[id==x]), as.numeric(methyl[id==x]), use = "pairwise.complete.obs")))
  }
})