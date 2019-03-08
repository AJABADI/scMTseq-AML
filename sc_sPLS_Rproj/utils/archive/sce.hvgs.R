sce.hvgs = function(sce,n.hvgs=1000){
  fit <- trendVar(sce, method="loess", use.spikes=FALSE)  ## fit a mean-dependent loess to variance
  var.decomp <- decomposeVar(sce, var.fit) ## decompose variance to technical and biological
  hvgs <- rownames(var.decomp[order(-var.decomp$bio)[1:n.hvgs],]) ## order them by biological variability
  return(hvgs)
}