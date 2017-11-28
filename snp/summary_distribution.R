#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)}

library(tools)

pi.fn <- args[1]
output_file <- paste(file_path_sans_ext(pi.fn), ".stat", sep="")

my_data <- read.csv(pi.fn, sep="\t", header=TRUE)
stats <- as.vector(summary(my_data$PI))
transposed <- t(stats)
write.table(transposed, file=output_file, sep="\t", row.names=FALSE, col.names=FALSE)