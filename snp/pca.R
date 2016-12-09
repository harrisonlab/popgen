#!/usr/bin/env Rscript

#Print a PCA plot calculated based on SNP matrix.
#First argument (required): input VCF file
#Output: a PDF file with the figure, suffix: _pca.pdf 

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)}

library("SNPRelate")
library("gdsfmt")
library("ggplot2")
library("ggrepel")
library(tools)

#The package offers also IBD and IBS calculations
vcf.fn <- args[1]
file_gds <- paste(file_path_sans_ext(vcf.fn), ".gds", sep="")
snpgdsVCF2GDS(vcf.fn, file_gds, method="biallelic.only")
snpgdsSummary(file_gds)
genofile <- snpgdsOpen(file_gds)

#PCA including haploid data
pca <- snpgdsPCA(genofile, autosome.only=FALSE)

#Percentage of variation explained by each axis
pc.percent <- pca$varprop*100
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1], # the first eigenvector
                  EV2 = pca$eigenvect[,2], # the second eigenvector
                  stringsAsFactors = FALSE)
##Assign populations identifiers to each individual
#tab$pop <- c("1", "1", "2", "1", "2", "2", "2", "1", "2", "2")

plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
variance <- (round(pc.percent,2))
pca_1 <- variance[1]
pca_2 <- variance[2]
xlabel <- paste("PC2 (", pca_2, "%)", sep="")
ylabel <- paste("PC1 (", pca_1, "%)", sep="")
pca_plot <- ggplot(tab, aes(EV2,EV1), fill="red") + geom_point(size=2) + geom_text_repel(aes(label=sample.id)) + xlab(xlabel) + ylab(ylabel)
#pca_plot <- ggplot(tab, aes(EV2,EV1), fill="red") + geom_point(colour=factor(tab$pop), size=2) + geom_text_repel(aes(label=sample.id)) + xlab(xlabel) + ylab(ylabel)
pca_plot2 <- pca_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
file_pca <- paste(file_path_sans_ext(vcf.fn), "_pca.pdf", sep="")
ggsave(file_pca, pca_plot2, dpi=300, height=5, width=5)
