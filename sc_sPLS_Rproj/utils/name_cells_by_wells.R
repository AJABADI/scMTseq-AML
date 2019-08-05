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
  underscores <- underscores+1 ## add one to adjust for addition of zerp
  substr_fun <- function(x,y, pos) substr(x, y[pos[1]]+1,y[pos[2]]-1)
  names <-  mapply(substr_fun,
                   x=char_vec,
                   y=pos_underscore, MoreArgs = list(pos=underscores))
  return(unname(names))
}

## in-line test

c01 <- c("A5_asas", "A9_shdghdg")
c12 <- c("x10_A61_hgh","ghsd_A5_hjhsd")

if(!all(name_cells_by_wells(c01, c(0,1))==c("A5","A9"), 
        name_cells_by_wells(c12, c(1,2))==c("A61", "A5"))){
  stop("name_cells_by_wells test failed")
}
