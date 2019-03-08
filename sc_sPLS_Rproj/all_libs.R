suppressMessages(library(mixOmics)) ## sPLS and PCA
suppressMessages(library(SingleCellExperiment)) ## sc assays
suppressMessages(library(data.table)) ## to use data.tables
suppressMessages(library(scran)) ## for sc analysis
suppressMessages(library(scater)) ## for QC
suppressMessages(library(Rtsne)) ## for tSNE
suppressMessages(library(ggplot2)); theme_set(theme_bw())
suppressMessages(library(ggrepel)) ## for nice ggplot labelling
suppressMessages(library(gridExtra)) ## to plot a grid pf pca plots
suppressMessages(library(grid)) ## to plot a grid pf pca plots
suppressMessages(library(umap)) ## for dimension reduction
suppressMessages(library(magrittr)) ## for pipe operator
suppressMessages(library(latex2exp)) ## for v_hat and r_bar labels
suppressMessages(library(stringr)) ## to use str_c for directories