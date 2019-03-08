#################### function to spot the well names using underscores b/w which the names can be found.
## INPUT:
    ## char_vec: character vector of unprocessed cell names
    ## underscores: the index of the underscores. e.g. for cell A101 with raw name 'abc_zx_A101_ghf_dy', it is c(2,3)
    ## as the name is b/w the 2nd and the 3rd underscores

name_cells_by_wells <- function(char_vec,
                               underscores = c(1,2)){ ## index of undersores after and before whose the well name exists
  pos_underscore <- gregexpr('_', char_vec)
  ## add the zeroth underscore at 0
  pos_underscore = lapply(pos_underscore, function(x) c(0,x))
  names <-  mapply(function(x,y) substr(x, y[1]+1,y[2]-1),
                   char_vec,
                   pos_underscore)
  return(unname(names))
}