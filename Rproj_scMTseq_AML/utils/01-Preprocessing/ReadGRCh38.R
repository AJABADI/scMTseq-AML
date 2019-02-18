exons = read.delim(file = '../data/scBS-seq/GRCh38/Exons.txt')
View(exons[1:10,])

introns = read.delim(file='data/scBS-seq/GRCh38/Introns.txt')
View(introns[1:10,])

enhancers = read.delim(file='data/scBS-seq/GRCh38/Enhancers_centred_2kb.txt')
View(enhancers[1:10,])

promoters = read.delim(file='data/scBS-seq/GRCh38/Promoter_mRNA_upstream_1500_500.txt')
View(promoters[1:10,])

CGI = read.delim(file='data/scBS-seq/GRCh38/CGI_SeqMonk.txt')
head(CGI[1:10,])

repMask =read.delim(file ='data/scBS-seq/GRCh38/repMask-LINE.L1.csv')
View(repMask[1:10,])
