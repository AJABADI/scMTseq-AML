## spls_out.RData

List of sPLS gene signature and their genomic ranges.


In all studies the link b/w expression and a context methylome (genebody, promoter, promoter cgi and promoter non cgi) is represented with a set of 50 genes from each dataset. It means that the variation of expression of RNAseq signature genes was most strongly linked to variation in BSseq signature genes for every context given.


The nested list gives: 

  - for every genomic context,
    -- for every sPLS component (1,2,3) (think of it as a PCA component, but summarising associations), and 3):
        --- for signature genes from each 'Omic, i.e. RNAseq or BSseq:
            1. The pairwise Pearson correlation violin plots and 
            2. A data.frame of the signature genes and their info in following format:

            -------------------------------------------------------------
            rowname              chr     start        end        strand
            ENSG00000041357       15    8540405      78552419       +
            -------------------------------------------------------------
            

There are some missing gene info in RNAseq signature which means the IDs did not match with the bed files which I can look into (It could be becase they're from chromosomes not included in bed files) .

Pointer: interesting ones that I looked into and the pairwise correlations had non-zero medians:

i) Promoter, componenet 1 , BSseq genes

* BSseq signature plots and gene info:
spls_out$prom$comp_1$BSseq (there's a $plots and a $signature slot): The methylation in the promoter regions from BSseq data whose methylation is linked to transcriptome heterogeneity tend to negatively correlate with the expression of their own genes.

ii) Genebody, component 1 (and possibly 2), BSseqs

* BSseq signature plots and gene info:
spls_out$genebody$comp_1$BSseq: The  methylation in the genebody regions from BSseq data whose methylation is linked to transcriptome heterogeneity tend to positively correlate with the expression of themsevles. For these genes, the correlation pattern b/w promoter methylation and expression is negative for those with promoters in CGI's and positive in the other ones.
