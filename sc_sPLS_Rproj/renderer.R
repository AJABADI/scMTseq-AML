## run this and all Rmd files are re-rendered and updated
## to update individual files, open them and knit

setwd('/Users/alabadi/Documents/_Projects/Cheatsheet/GitHub')
rmds <- list.files(pattern = ".Rmd", recursive = TRUE, full.names = TRUE)
## exclude index and template
rmds %<>% .[!basename(rmds) %in% c('foo-cheatsheet.Rmd', 'index.Rmd' )] ## index should be rendered at the end

exclude <- c() # c("./01-computer/R/R-cheatsheet.Rmd")
for (i in rmds){
  ## exclude if desired
  if(!i %in% exclude){
    html_name <- i %>% basename() %>% str_replace(., 'Rmd', 'html')
    rmarkdown::render(input = i, output_format = 'html_document', params = NULL, output_file = html_name, output_dir = 'docs' )
  }
}

rmarkdown::render(input = "index.Rmd", output_format = 'html_document', params = NULL, output_file = "index.html", output_dir = 'docs' )



# 
# rmarkdown::render(input, output_format = NULL, output_file = NULL, output_dir = NULL,
#                   output_options = NULL, intermediates_dir = NULL,
#                   knit_root_dir = NULL,
#                   runtime = c("auto", "static", "shiny", "shiny_prerendered"),
#                   clean = TRUE, params = NULL, knit_meta = NULL, envir = parent.frame(),
#                   run_pandoc = TRUE, quiet = FALSE, encoding = getOption("encoding"))
