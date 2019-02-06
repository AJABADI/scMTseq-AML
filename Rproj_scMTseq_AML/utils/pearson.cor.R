## function to calculate pearson correlation b/w expressio and methylation at each genomic context
source('utils/calc.cv.R')
pearson.cor = function(rna.sce=rna.sce, met.sce=met.gene.sce.filt_NAs, n.hvgs=5000, calc.cv = calc.cv){
  ## get the counts
  expr = logcounts(rna.sce)
  met = counts(met.sce)
  ## intersect genes
  intersect = Reduce(intersect , list(rownames(expr), rownames(met)))
  expr = expr[intersect,]
  met=met[intersect,]
  if (length(intersect)<n.hvgs){
    n.hvgs = length(intersect)
    cat("Number of hvgs overridden by common genes : ", length(intersect))
  }
  ## get HVGs in methylation data
  hvgs.met = calc.cv(met, n.hvgs = n.hvgs)$hvgs
  ## keep only these HVGs in both datasets
  expr.hvg = expr[hvgs.met,]
  met.hvg=met[hvgs.met,]
  pearson = 0
  ## calculate the pearson corr for pairwise complete observations
  for (i in 1:dim(met.hvg)[1]){
    pearson[i] = cor(expr.hvg[i,], met.hvg[i,], use = "pairwise.complete.obs")
  }
  names(pearson) = rownames(met.hvg)
  
  out=list()
  
  out$pearson = pearson
  ## calculate the variation of methylation and expression for each gene
  all.genes = data.frame(expr = calc.cv(expr)$df$vars,
                       met = calc.cv(met)$df$vars)
  rownames(all.genes) = rownames(expr)
  ## calculate the variation of methylation and expression for hvgs
  hv.genes = data.frame(expr = calc.cv(expr.hvg)$df$vars,
                       met = calc.cv(met.hvg)$df$vars)
  rownames(hv.genes) = rownames(expr.hvg)
  
  vars = list('hvgs' = hv.genes, 'all' = all.genes)
  out$vars = vars
  return(out)
}