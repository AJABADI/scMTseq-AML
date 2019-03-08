#################### GLOBAL PARAMETERS
params_global = list(
  # data_dir = '/Volumes/0414861919/Mac_Home/Documents/Projects/sc_sPLS/sc_sPLS',
  data_dir = '../',
  out_dir= '../output',
  RData = 'RData',
  save.RData = !FALSE,
  load.RData= TRUE, ## for time-taking ones only
  eval=TRUE, ## FALSE: load session data and just render
  cache=TRUE,
  echo=!TRUE
)
## check
stopifnot(all(
  with(params_global,
       all(dir.exists(out_dir),
           dir.exists(RData),
           any(save.RData&eval,!save.RData)))
))
#################### 01-scRNAseq PARAMETERS
params_01= list(
  ## raw RNA data
  rna_file= file.path(params_global$data_dir, 'data/scRNAseq/raw/Kasumi1_scRNAseq_RawCounts.txt'),
                      ## gene filtering
                      min.gene.libsize = log10(4), ## minimum log expression
                      max.cell.dropout = 0.95, ## maximum level of cell dropout
                      rna.78cells.hvgs = 2000 ## hvgs
  )
  ## check
  stopifnot(all(
    with(params_01,
         all(file.exists(rna_file)))
  ))
  #################### 02-scBSseq PARAMETERS
  params_02= list(
    ## raw single bp met calls folder
    met_1bp_in= file.path(params_global$data_dir, 'data/scBSseq/met/raw/bismark'),
    ## processed single bp met calls file name
    met_1bp_out= file.path(params_global$data_dir, 'data/scBSseq/met/parsed/met_data.tsv'),
    ## whether to use parsed data or (re-)parse
    use.parsed= TRUE,
    ## 3kb windows data folder
    met_3kb= file.path(params_global$data_dir, 'data/scBSseq/met/parsed/methylVariance-w3K-s1K5.csv'),
    ## .bed files for summarissing calls on genomic features
    bedfiles= file.path(params_global$data_dir, 'data/scBSseq/filt'),
    cov.cutoff = 0.5, ## coverage cut-off filtering for genes
    hvgs_pca = 2000 ## hvgs to use for pca
  )
  
  ## check
  stopifnot(all(
    with(params_02,
         all(
           file.exists(c(met_1bp_in, met_3kb)),
           dir.exists(c(substr(met_1bp_out, 0, gregexpr("\\/[^\\/]*$",met_1bp_out)[[1]][1]), bedfiles)),
           !use.parsed |(use.parsed&file.exists(met_1bp_out))
         ))
  ))
  #################### 02-scBSseq PARAMETERS
  params_03= list(
    scale = FALSE, ## scale datasets?
    fig_asp = 0.7 ## spls plot aspect ratios
  )
  ####################
  params = c(params_global, params_01, params_02, params_03)
  