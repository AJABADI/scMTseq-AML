## seurat clustering
library(Seurat)
seurat <- CreateSeuratObject(counts = counts(sce_rna), project = "pbmc3k", min.cells = 3, min.features = 200)
nbt=setup(nbt,project="NBT",min.cells = 3,names.field = 2,names.delim = "_",min.genes = 1000,is.expr=1,)
Seurat::setu