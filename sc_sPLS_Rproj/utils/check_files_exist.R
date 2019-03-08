## function to check a list of entry files exist and if not, produce and warning mentionig the files that do not exist

check_files_exist = function(file_list = io$files){
  non_exist = lapply(file_list, file.exists) %>% .[!unlist(.)] %>% names()
  if(length(non_exist)) warning(paste0('file does not exist: ', non_exist))
}