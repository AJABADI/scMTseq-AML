block_spls_plot = function(block = block_sPLS, ## which elements in the list?
                     pars = par_plot,
                     x=c(2,3),
                     y=1
){
  
  shapes = par_plot$shape[c(x,y)]
  cols = par_plot$col[c(x,y)] %>%  as.character()
  labs = as.character(par_plot$label)[c(x,y)]
  
    ## correlation circle plot
    plotVar(block, legend.title = "", legend = TRUE ,pch=shapes,
            col= cols, title ="")
    ## sample plot
    plotIndiv(block, subtitle = LETTERS[1:3], size.subtitle = 12)
    ## arrow plot
    plotArrow(block, ind.names = TRUE, X.label = "", Y.label = "", title = "Block sPLS")
  
}