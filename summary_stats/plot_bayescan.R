#!/usr/bin/env Rscript

source('/home/sobczm/bin/bayescan2.1/R_functions/plot_R.r')
library(tools)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file)", call.=FALSE)}
input <- args[1]

out <- plot_bayescan(input,0,FDR=0.05)
out_table <- paste(file_path_sans_ext(input), "_outliers.txt", sep="")
write.table(out$outliers, file=out_table, sep ="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

y <- read.table(input)

#Plot the distribution of FST values
fst <- y[,5]
d <- density(fst)
out_fst <- paste(file_path_sans_ext(input), "_fst.pdf", sep="")
pdf(out_fst, width=11, height=8)
plot(d, main="Fst")
polygon(d, col="purple", border="black") 
dev.off()

#Plot the distribution of alpha coefficients
alpha <- y[,4]
d <- density(alpha)
out_alpha <- paste(file_path_sans_ext(input), "_alpha.pdf", sep="")
pdf(out_alpha, width=11, height=8)
plot(d, main="alpha")
polygon(d, col="blue", border="black") 
dev.off()

#Plot the distribution of qvalues
qval <- y[,3]
d <- density(qval)
out_qval <- paste(file_path_sans_ext(input), "_qval.pdf", sep="")
pdf(out_qval, width=11, height=8)
plot(d, main="qval")
polygon(d, col="green", border="black") 
dev.off()
