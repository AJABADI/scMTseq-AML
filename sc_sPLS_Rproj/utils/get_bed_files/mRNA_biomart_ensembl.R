
# This script allows us to extract the following information from an ENSEMBL database:
# - Genes
# - Transcripts
# - Exons
# - UTRs (all and 5' and 3' separately)
# - Coding sequences (coding exons)

library(stringr)
library(dplyr)
library(biomaRt)
listMarts(archive=F)
source("/Users/ricard/data/ensembl/utils.R")

# outdir <- "/Users/ricard/data/ensembl/human/v87/BioMart/mRNA"; dir.create(outdir, recursive=T)
outdir <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA"; dir.create(outdir, recursive=T)

# Define database and dataset
version = "87"
# version = "75"
database <- "ENSEMBL_MART_ENSEMBL"
specie = "Mmusculus"
# host="may2012.archive.ensembl.org"
host="dec2016.archive.ensembl.org"

if (specie == "Hsapiens") {
  allchr = c(1:22,"X","Y","MT")
  dataset <- "hsapiens_gene_ensembl"
  gene_attributes <- c("chromosome_name", "start_position", "end_position","strand",
                       "ensembl_gene_id", "description", "hgnc_symbol")
}
if (specie == "Mmusculus") {
  allchr = c(1:19,"X","Y","MT")
  dataset <- "mmusculus_gene_ensembl"
  gene_attributes <- c("chromosome_name", "start_position", "end_position","strand",
                       "ensembl_gene_id", "description", "mgi_symbol")
}

ensembl <- useMart(biomart=database, dataset=dataset, verbose=T, host=host)

## Retrieve protein-coding gene information ##
genes = getBM(attributes=gene_attributes, mart = ensembl,
              # filters=c("chromosome_name","biotype"), values=list(allchr,"protein_coding")) %>% tbl_df
              filters=c("chromosome_name"), values=list(allchr)) %>% tbl_df
colnames(genes) <- c("chr","start","end","strand","ens_id","description","symbol")
genes$chr <- convert_chr_format(genes$chr, to="long")
genes$strand <- c("-","+")[factor(genes$strand)]
genes <- genes[with(genes, order(chr,start,end)), ]
genes <- remove(genes, riken=T, no_symbol=T, dna_segments=T, empty_description=T, duplicated=T, pseudogenes=T, predicted=T, expressed_sequence=T, olfactory=T)
write.table(genes,str_c(outdir,"/",specie,"_genes_BioMart.",version,".txt",sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)


## Retrieve transcript information ##
transcript_attributes <- c("chromosome_name", "start_position", "end_position","strand","ensembl_gene_id",
                           "transcript_start", "transcript_end", "ensembl_transcript_id","transcript_biotype")
transcripts = getBM(attributes=transcript_attributes, mart=ensembl,
                    filters=c("chromosome_name","biotype"), values=list(allchr,"protein_coding")) %>% tbl_df
colnames(transcripts) <- c("chr","gene_start","gene_end","strand","ens_gene_id","start","end","ens_transcript_id","transcript_biotype")
transcripts$chr <- convert_chr_format(transcripts$chr, to="long")
transcripts$strand <- c("-","+")[factor(transcripts$strand)]
transcripts <- transcripts[with(transcripts, order(chr,start,end)), ]
transcripts <- filter(transcripts, ens_gene_id %in% genes$ens_id)
transcripts <- filter(transcripts, transcript_biotype=="protein_coding") %>% dplyr::select(-transcript_biotype)
stopifnot(all(transcripts$ens_gene_id %in% genes$ens_id))
write.table(transcripts,str_c(outdir,"/",specie,"_transcripts_BioMart.",version,".txt",sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)


## Retrieve exon information ##
exon_attributes <- c("chromosome_name", "start_position", "end_position", "strand", "ensembl_gene_id",
                     "transcript_start", "transcript_end", "ensembl_transcript_id","transcript_biotype",
                     'ensembl_exon_id', "exon_chrom_start", "exon_chrom_end")
exons <- getBM(attributes=exon_attributes, mart=ensembl,filters=c("chromosome_name","biotype"), values=list(allchr,"protein_coding")) %>% tbl_df
colnames(exons) <- c("chr","gene_start","gene_end","strand","ens_gene_id","transcript_start","transcript_end","ens_transcript_id","transcript_biotype","ens_id","start","end")
exons$chr <- convert_chr_format(exons$chr, to="long")
exons$strand <- c("-","+")[factor(exons$strand)]
exons <- exons[with(exons, order(chr,start,end)), ]
exons <- filter(exons, transcript_biotype=="protein_coding") %>% dplyr::select(-transcript_biotype)
exons <- filter(exons, ens_gene_id %in% genes$ens_id)
write.table(exons,str_c(outdir,"/",specie,"_exons_BioMart.",version,".txt",sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)


## Retrieve coding sequence information ##
cds_attributes <- c("chromosome_name", "start_position", "end_position", "strand", "ensembl_gene_id",
                    "transcript_start", "transcript_end", "ensembl_transcript_id","transcript_biotype",
                    'ensembl_exon_id',"exon_chrom_start", "exon_chrom_end","genomic_coding_start","genomic_coding_end")
cds <- getBM(attributes=cds_attributes, mart=ensembl, filters=c("chromosome_name","biotype"), values=list(allchr,"protein_coding")) %>% tbl_df

colnames(cds) <- c("chr","gene_start","gene_end","strand","ens_gene_id","transcript_start","transcript_end","ens_transcript_id","transcript_biotype","ens_exon_id","exon_start","exon_end","cds_start","cds_end")
cds$chr <- convert_chr_format(cds$chr, to="long")
cds$strand <- c("-","+")[factor(cds$strand)]
cds <- cds[with(cds, order(chr,gene_start,gene_end)), ]
cds <- filter(cds, transcript_biotype=="protein_coding") %>% dplyr::select(-transcript_biotype)
cds <- cds %>% filter(!(is.na(cds_start) & is.na(cds_end)))
cds <- filter(cds, ens_gene_id %in% genes$ens_id)
stopifnot(all(cds$ens_gene_id %in% genes$ens_id))
write.table(cds,str_c(outdir,"/",specie,"_cds_BioMart.",version,".txt",sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)


## Extract UTRs ##
utr_attributes <- c("chromosome_name", "start_position", "end_position","strand","ensembl_gene_id",
                    "transcript_start", "transcript_end", "ensembl_transcript_id","transcript_biotype",
                    "5_utr_start", "5_utr_end", "3_utr_start", "3_utr_end")
utr <- getBM(attributes=utr_attributes, mart=ensembl, filters=c("chromosome_name","biotype"), values=list(allchr,"protein_coding")) %>% tbl_df
colnames(utr) <- c("chr","gene_start","gene_end","strand","ens_gene_id","transcript_start","transcript_end","ens_transcript_id","transcript_biotype","five_utr_start","five_utr_end","three_utr_start","three_utr_end")
utr$chr <- convert_chr_format(utr$chr, to="long")
utr$strand <- c("-","+")[factor(utr$strand)]
utr <- utr[with(utr, order(chr,five_utr_start,five_utr_end)), ]
utr <- filter(utr, transcript_biotype=="protein_coding") %>% dplyr::select(-transcript_biotype)
utr <- utr %>% filter(! (is.na(three_utr_start) & is.na(five_utr_start) & is.na(five_utr_end) & is.na(three_utr_end) ))
utr <- filter(utr, ens_gene_id %in% transcripts$ens_gene_id)
utr <- filter(utr, ens_transcript_id %in% transcripts$ens_transcript_id)
stopifnot(all(utr$ens_gene_id %in% genes$ens_id))
# write.table(five_utr,str_c(outdir,"/Hsapiens_utr_BioMart.",version,".txt",sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

# 5'UTR
five_utr <- utr %>% filter(!(is.na(five_utr_start) & is.na(five_utr_end))) %>% dplyr::select(-contains(c("three_utr")))
colnames(five_utr) <- c("chr","gene_start","gene_end","strand","ens_gene_id","transcript_start","transcript_end","ens_transcript_id","start", "end")
write.table(five_utr,str_c(outdir,"/",specie,"_5utr_BioMart.",version,".txt",sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

# 3'UTR
three_utr <- utr %>% filter(!(is.na(three_utr_start) & is.na(three_utr_end))) %>% dplyr::select(-contains(c("five_utr")))
colnames(three_utr) <- c("chr","gene_start","gene_end","strand","ens_gene_id","transcript_start","transcript_end","ens_transcript_id","start", "end")
write.table(three_utr,str_c(outdir,"/",specie,"_3utr_BioMart.",version,".txt",sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)


