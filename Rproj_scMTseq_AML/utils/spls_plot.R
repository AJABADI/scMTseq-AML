spls_plot = function(spls_list = spls_all, ## which elements in the list?
                     pars = par_plot,
                     hide = NULL,
                     hide_in_plotVar =FALSE, ## whether to exclude hide cells in plotVar as well
                     plots = c("plotVar", "plotIndiv", "plotArrow")
){
  
  Y_lab = pars[1,]$label %>% as.character()
  Y_shape = pars[1,]$shape
  Y_col = pars[1,]$col %>% as.character()
  
  for (i in rownames(pars[-1,])){
    spls_obj = spls_list[[i]]
    par = pars[i,]
    X_lab = paste0(par$label, " Methylome")
    col_XY = c(as.character(par$col), Y_col)
    pch_XY = c(par$shape, Y_shape)

    paste0("sPLS plots for ", Y_lab, " & ", X_lab, " :")
    ## correlation circle plot
    if((!hide_in_plotVar | is.null(hide)) & "plotVar" %in% plots){
      plotVar(spls_obj, legend.title = "", legend = c(X_lab, Y_lab) ,pch=pch_XY,
              col= col_XY, title ="" )
    }
    
    if(!is.null(hide)){
      spls_obj$names$sample %<>% .[!. %in% hide]
      spls_obj$variates$X %<>% .[!rownames(.) %in% hide,]
      spls_obj$variates$Y %<>% .[!rownames(.) %in% hide,]
      spls_obj$X %<>% .[!rownames(.) %in% hide,]
      spls_obj$Y %<>% .[!rownames(.) %in% hide,]
      
      if(hide_in_plotVar & "plotVar" %in% plots){
        plotVar(spls_obj, legend.title = "", legend = c(X_lab, Y_lab) ,pch=pch_XY,
                col= col_XY, title ="" )
      }
    }
    if("plotIndiv" %in% plots){
      ## sample plot
      plotIndiv(spls_obj, subtitle = c(X_lab, Y_lab), size.subtitle = 12)
    }
    if("plotArrow" %in% plots){
      ## arrow plot
      plotArrow(spls_obj, ind.names = TRUE, X.label = X_lab, Y.label = Y_lab)
    }
  }
  
}