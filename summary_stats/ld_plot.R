#!/usr/bin/env Rscript

#Plots r^2 (measure of linkage disequilibrium) as a heatmap.
########################################################################
#Input: VCFtools output from the command --hap-r2.
########################################################################
#Outputs: #Print R2 LD plot for each contig, in the PNG (low resolution)
#and PDF (high resolution) formats.
#########################################################################
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)}

library(gg.ldplot)
table <- args[1]
vcf_ld <- read.table(table, header = TRUE, stringsAsFactors = F)
contig_list <- split(vcf_ld, vcf_ld$CHR)


for (x in 1:length(contig_list))
{
filename <- paste(contig_list[[x]][1,1], "_ld_plot.pdf", sep="")
filename2 <- paste(contig_list[[x]][1,1], "_ld_plot.png", sep="")
pdf(filename, width=7, height=5)
print(plot_LDmatrix(contig_list[[x]]))
dev.off()
png(filename2, width=7, height=5)
print(plot_LDmatrix(contig_list[[x]]))
dev.off()
}