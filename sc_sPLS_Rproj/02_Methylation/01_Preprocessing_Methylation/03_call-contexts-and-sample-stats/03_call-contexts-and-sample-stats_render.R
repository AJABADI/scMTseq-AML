## for HPC
.libPaths('~/R_libs/')
## check Rmd parameters before full run
R_file_wd <- file.path(getwd(),'02_Methylation/01_Preprocessing_Methylation/03_call-contexts-and-sample-stats')
rmd_file <- '03_call-contexts-and-sample-stats.Rmd'
test_run <- FALSE

## if not on hpc
if(grepl('Al',Sys.info()['nodename'])){
  setwd(R_file_wd) ##CHECK
}

out_dir=NULL
docs <- '../../../docs'
if(dir.exists(docs)){
  out_dir <- docs
}

rmarkdown::render(input = rmd_file,  ##CHECK the rmd file
                  # output_format = 'html_document', NULL==default
                  # output_file = "index.html", ## name of output html - NULL==default
                  # params = list(cache=FALSE, echo=TRUE, 
                  #               chrs=c(1:22,"MT"), ##CHECK chromosomes to use for genome c stats
                  #               genome_stats_rds='genome_c_stats.Rds'), ##CHECK
                  output_dir = out_dir ## directory of output html
)
