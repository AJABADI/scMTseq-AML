circosPlot_spls <- function(spls.obj,  na.rm=TRUE, cutoff=0.7, ...){
  ## circosPlot uses X,Y,indY,loadings, varites, ncomp
  Y=rep_len(1:2, dim(spls.obj$Y)[1])
  X=list(X=spls.obj$X, Y=spls.obj$Y)
  ncomp=spls.obj$ncomp
  list.keepX = list(X = rep(10, ncomp), Y = rep(10,ncomp))
  blocksplsda <- block.splsda(X = X, Y = Y,
                              ncomp = ncomp, keepX = list.keepX)
    
  blocksplsda$X$X <- spls.obj$X
  blocksplsda$X$Y <- spls.obj$Y
  
  blocksplsda$variates$X <- spls.obj$variates$X
  blocksplsda$variates$Y <- spls.obj$variates$Y
  
  blocksplsda$loadings$X <- spls.obj$loadings$X
  blocksplsda$loadings$Y <- spls.obj$loadings$Y
  
  if(na.rm){
    data <- blocksplsda$X[[1]]
    for(i in 1:ncol(data)){
      data[is.na(data[,i]), i] <- mean(data[,i], na.rm = TRUE)
    }
    data -> blocksplsda$X[[1]]
    rm(data)
  }
  circosPlot(blocksplsda,cutoff=cutoff, ...)
}
