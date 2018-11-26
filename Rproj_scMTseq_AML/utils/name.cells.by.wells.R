#################### function to spot the well names using underscores
## INPUT:
    ## char.vec: character vector of unprocessed cell names
    ## pos.in.string: whether the cell name is at first or middle of string

name.cells.by.wells = function(char.vec,
                               pos.in.string = 'middle'){ ## position = c('first'. 'middle')
  pos.underscore = gregexpr('_', char.vec)
  
  for (i in 1:length(char.vec)){
    
    if(pos.underscore[[i]][1]>1){
      if(pos.in.string =='first'){ ## if well name at the beginning of string
        char2 = pos.underscore[[i]][1] - 1
        char.vec[i] = substr(char.vec[i],0, char2)
      } else { ## if well name in the middle of the string
        char1 = pos.underscore[[i]][1] + 1
        char2 = pos.underscore[[i]][2] - 1
        char.vec[i] = substr(char.vec[i],char1, char2)
      }
    }
  }
  return(char.vec)
}