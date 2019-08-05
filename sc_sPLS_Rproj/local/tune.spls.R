suppressMessages({
library(data.table)
library(mixOmics)
})

opts <- list()
opts$min_cov_full <- 25
opts$min_cov_filt <- 22
opts$filt.cells <- c('H8', 'E8', 'H6')
keeps <- c(50,50)

# regression full set ----
list_spls$regression$fullset <- lapply(annos, FUN = function(x) 
  spls_fun(met_dt = met_exp, anno = x, min.cov = opts$min_cov_full,
           mode='regression', filt.cells=NULL, Y=transcriptome, keepX=keeps, keepY=keeps))

# regression filtered ----
list_spls$regression$filt <- lapply(annos, FUN = function(x)
  spls_fun(met_dt = met_exp, anno = x, min.cov = opts$min_cov_filt,
           mode='regression', filt.cells=opts$filt.cells, Y=transcriptome, keepX=keeps, keepY=keeps))

# canonical full set ----
list_spls$regression$fullset <- lapply(annos, FUN = function(x)
  spls_fun(met_dt = met_exp, anno = x, min.cov = opts$min_cov_full,
           mode='canonical', filt.cells=NULL, Y=transcriptome, keepX=keeps, keepY=keeps))

# canonical filtered ----
list_spls$regression$filt <- lapply(annos, FUN = function(x)
  spls_fun(met_dt = met_exp, anno = x, min.cov = opts$min_cov_filt,
           mode='canonical', filt.cells=opts$filt.cells, Y=transcriptome, keepX=keeps, keepY=keeps))