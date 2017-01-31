#!/usr/bin/env Rscript

#Plot the output of linkage diseqilibrium (D, D'. r2) calculations carried out by
#by vcftools (script fus_linkage_disequilibrium.sh)
#First argument (required): input file with results table output by vcftools
#Output: Scatterplots - r2 (file suffix _r2.pdf) and D' (file suffix _d.pdf) versus physical distance
#between the markers; Histogram of D' values (file suffix _d_hist.pdf )

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)}

library(ggplot2)
library(tools)
table <- args[1]
input <- as.data.frame(read.table(table,header=TRUE))
distance <- input$POS2 - input$POS1
#Plot R^2 against distance
r2_d <- paste(file_path_sans_ext(table), "_r2.pdf", sep="")
pdf(r2_d, width=11, height=8)
qplot(distance, input$R.2)
dev.off()
#Plot D' against distance
d_d <- paste(file_path_sans_ext(table), "_d.pdf", sep="")
pdf(d_d, width=11, height=8)
qplot(distance, input$Dprime)
dev.off()
#And histogram
hist_d <- paste(file_path_sans_ext(table), "_d_hist.pdf", sep="")
pdf(hist_d, width=11, height=8)
qplot(input$Dprime, geom="histogram")
dev.off()
