## function to perform (block) sPLS on group of methylation sces against
## transcription sce given proportion/number of top genes required in each sce

spls_wrapper = function(df = par_plot, x = 3, y = 1, ncomp =2,
                        top_x = NULL, top_y = NULL, cov = NULL, exclude=NULL,
                        near_zer_var =TRUE, ## if near zero variance is expected
                        keepX=c(50,50), keepY=c(50,50), block=FALSE, mode="canonical"){
  
  ## Y block sce
  y_sce = get(rownames(df[y,]))
  ## cell filtering
  cells = colnames(y_sce)
  if(!is.null(exclude)){ 
    cells = cells[!cells %in% exclude]
    y_sce = y_sce[,cells]
  }
  
  ## Y gene filtering by variation
  if(!is.null(top_y)){
    top_y = ifelse(top_y>1,min(dim(y_sce)[1], top_y),
                     floor(dim(y_sce)[1]*top_y))
    y_sce = y_sce[1:top_y,]
  }
    
  Y = t(logcounts(y_sce))
  keep_X = keepX
  
  if(block){
    if(length(x)<=1) stop("x must be a vector for block spls")
    X = list()
    keep_X = list()
    for (i in x ){
      sce_name = rownames(df[i,])
      x_sce = get(sce_name)
      
      ## gene filtering by variation
      if(!is.null(top_x)){
        top_x = ifelse(top_x>1,min(dim(x_sce)[1], top_x),
                         floor(dim(x_sce)[1]*top_x))
        x_sce = x_sce[1:top_x,]
      }
      
      ## gene filtering by coverage
      if(!is.null(cov)){
        x_sce = x_sce[rowData(x_sce)$cov>=cov,]
        message(paste0("Number of genes remained in methylation: ", dim(x_sce)[1]))
      }
      
      
      X[[sce_name]] = t(assay(x_sce, "rates"))
      keep_X[[sce_name]] = keepX
      if(!is.null(exclude)){ 
        X[[sce_name]] %<>% .[cells,]
      }
    }

    out = block.spls(X=X, Y=Y,ncomp=ncomp ,keepX=keep_X, 
                     keepY=keepY, mode=mode, near.zero.var = near_zer_var)
  } else {
    out = list()
    for (i in x){
      sce_name = rownames(df[i,])
      x_sce = get(sce_name)
      if(!is.null(exclude)){ 
        x_sce = x_sce[,cells]
      }
      X = t(assay(x_sce, "rates"))
      
      ## gene filtering by variation
      if(!is.null(top_x)){
        top_x = ifelse(top_x>1,min(dim(x_sce)[1], top_x),
                       min(dim(x_sce)[1], floor(dim(x_sce)[1]*top_x)))
        x_sce = x_sce[1:top_x,]
      }
      
      ## gene filtering by coverage
      if(!is.null(cov)){
        x_sce = x_sce[rowData(x_sce)$cov>=cov,]
        message(paste0("Number of genes remained in methylation: ", dim(x_sce)[1]))
      }
      
      X = t(assay(x_sce,"rates"))
      
      out[[sce_name]] = spls(X=X, Y=Y,ncomp=ncomp ,keepX=keep_X, 
                 keepY=keepY, mode=mode, near.zero.var = near_zer_var)
    }
    
  }

  return(out)
}