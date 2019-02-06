params = list(
  ## raw RNA data
  rna_file= '../data/scRNA-seq/raw/Kasumi1_scRNAseq_RawCounts.txt',
  ## raw single bp met calls folder
  met_1bp_in= '../data/scBS-seq/met/1bp/raw',
  ## processed single bp met calls file name
  met_1bp_out= '../data/scBS-seq/met/1bp/parsed/met_data.tsv',
  ## whether to use parsed data or (re-)parse
  use.parsed= F,
  ## 3kb windows data folder
  met_3kb= '../data/scBS-seq/met/3kbp/raw/methylVariance-w3K-s1K5.csv',
  ## .bed files for summarissing calls on genomic features
  bedfiles= '../data/scBS-seq/filt',
  RData= '../data/RData',
  save.RData = FALSE,
  load.RData= TRUE, ## for time-taking ones only
  eval=TRUE, ## FALSE: load session data and just render
  cache=TRUE,
  cov.cutoff = 0.5,
  echo=TRUE,
  
  ## gene filtering
  min.gene.libsize = log10(4), ## minimum log expression
  max.cell.dropout = 0.95, ## maximum level of cell dropout
  
  ## hvgs
  rna.78cells.hvgs = 2000,
    ## methylation
  hvgs_pca = 2000
)

with(params, stopifnot(all(max.cell.dropout >=0 & max.cell.dropout <=1,
                           cov.cutoff>=0&cov.cutoff<=1
                           )))

# with(params, stopifnot())


io = params
