## function to create sce from a given annotation in a data.table and
## also compute the coverage of cells by each gene in rowData
## DT to be of this form:
#      sample       id       anno     rate  weight
# 1:     A5 ENSG00000000003 genebody   29      7
# 2:     A5 ENSG00000000457 genebody  100     36
# 3:     A5 ENSG00000000460 genebody   91     43
# 4:     A5 ENSG00000000971 genebody   91     34
# 5:     A5 ENSG00000001036 genebody  100      1
# 6:     A5 ENSG00000001084 genebody   38    197


met_sce = function(met=met_all, ## the data.table of the above format
                   annot="genebody", ## one of unique(met$anoo)
                   cov = 0){
  ## entry check
  stopifnot(all(inherits(met, "data.table"),
                inherits(annot, "character"),
                # is.null(tot_c) |inherits(tot_c, "character"),
                cov>=0&cov<=1))
  run_time <- system.time({ ## record run time
  ## subset by annotation
  met_anno <- met[anno==annot][,anno:=NULL]
  
  
  ## subset by c context
  met_anno_cpg <- met_anno[c_context=="CpG"][,c_context:=NULL]
  met_anno_cph <- met_anno[c_context=="CpH"][,c_context:=NULL]
  
  ## weighted variance
  vhat = met_anno[,lapply(.SD,function(x)(weightedVar(as.numeric(x),w=weight))),by=id,.SDcols='rate']
  colnames(vhat)[2] = "vhat"
  ## weighted rate
  rbar = met_anno[,lapply(.SD,weighted.mean,w=weight),by=id,.SDcols='rate']
  colnames(rbar)[2] = "rbar"
  ## merge and make df
  wtd_arith = merge(rbar, vhat) %>% data.frame(row.names = "id")
  ## create weight and rate data.frames from the data.table
  rates = dcast(met[anno==annot], id~sample, value.var = 'rate') %>% data.frame(row.names = "id")
  weights = dcast(met[anno==annot], id~sample, value.var = 'weight')%>% data.frame(row.names = "id")
  ## modify this when you analyse windows - total number of possible methylation calls
  # if(!is.null(tot_c)){
  #   c_content = dcast(met[anno==annot], id~sample, value.var = tot_c)
  # }
  ## check that rownames are identical (just in case)
  stopifnot(all(identical(rownames(rates), rownames(wtd_arith)),
            identical(rownames(weights), rownames(wtd_arith) )))
  ## create sce
  sce = SingleCellExperiment(assays =SimpleList(rates = rates,
                                                 weights = weights), rowData = wtd_arith)
  ## order by weighted variance (and)
  sce %<>% .[order(-rowData(.)$vhat, -rowData(.)$rbar),]
  ## add dimension metadata
  metadata(sce)$gene_stats = data.frame()
  ## record original sce dimensions
  metadata(sce)$gene_stats["pre_filter",annot] = dim(sce)[1]
  
  ## add cell coverage proportion as rowData (proportion of cells in which the gene has had measurements)
  rowData(sce)$cov = rowSums(!is.na(assay(sce, "rates")))/ncol(sce)
  
  ## filter by coverage
  sce %<>% .[rowData(.)$cov>=cov,]
  ## add metadata for post-filtering
  metadata(sce)$gene_stats[paste0("cov_",round(cov,2)),annot] = dim(sce)[1]
  })
  ## ad run time to metadata
  metadata(sce)$run_time_summarise= unname(run_time[1])
  return(sce)
}
  

