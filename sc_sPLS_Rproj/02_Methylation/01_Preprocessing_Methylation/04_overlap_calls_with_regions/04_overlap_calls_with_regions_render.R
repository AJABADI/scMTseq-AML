## for HPC
.libPaths('~/R_libs/')
# rsync -avz /Users/alabadi/Documents/_Projects/sc_sPLS/sc_sPLS_Rproj/02_Methylation/01_Preprocessing_Methylation/04_overlap_calls_with_regions/04_overlap_calls_with_regions_render.R ajabadi@spartan.hpc.unimelb.edu.au:/data/cephfs/punim0613/AL/sc_sPLS/sc_sPLS_Rproj/02_Methylation/01_Preprocessing_Methylation/04_overlap_calls_with_regions/04_overlap_calls_with_regions_render.R
## check Rmd parameters before full run
R_file_wd <- file.path(getwd(),'02_Methylation/01_Preprocessing_Methylation/04_overlap_calls_with_regions')
rmd_file <- '04_overlap_calls_with_regions.Rmd'
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
