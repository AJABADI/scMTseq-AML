## if not on hpc
if(grepl('Al',Sys.info()['nodename'])){
  ## must be open in a project
  setwd(file.path(getwd(),'02_Methylation/01_Preprocessing_Methylation/01_genome_cytosine_stats')) ##INPUT
}
file_name <- '01_genome_cytosine_stats' ## Rmd file name ##INPUT
out_dir=NULL
docs <- '../../../docs'
if(dir.exists(docs)){
  out_dir <- docs
}

.libPaths('~/R_libs/')
rmarkdown::render(input = '01_genome_cytosine_stats.Rmd',  ##INPUT the rmd file
                  # output_format = 'html_document', NULL==default
                  # output_file = "index.html", ## name of output html - NULL==default
                  # params = list(cache=FALSE, echo=TRUE, 
                  #               chrs=c(1:22,"MT"), ##INPUT chromosomes to use for genome c stats
                  #               genome_stats_rds='genome_c_stats.Rds'), ##INPUT
                  output_dir = out_dir ## directory of output html
                  )
