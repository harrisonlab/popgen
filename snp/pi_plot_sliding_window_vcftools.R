#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)}
library(tools)
library(ggplot2)

pi.fn <- args[1]
my_data <- read.table(pi.fn, sep="\t", header=TRUE)
output_file <- paste(file_path_sans_ext(pi.fn), ".pdf", sep="")

ids <- max(my_data$BIN_END)
start <- min(my_data$BIN_START)
xaxis <- seq(from = start, to = ids, by = 200000)

slide_plot <- ggplot(my_data, aes(x=my_data$BIN_START, y=my_data$PI)) + geom_smooth(colour="black", fill="darkviolet") + scale_x_continuous(breaks = pretty(my_data$BIN_START, n = 20)) + xlab("Chrom coordinate") + ylab("Pi") + expand_limits( y = 0)

ggsave(output_file, slide_plot, dpi=300, height=6, width=20)