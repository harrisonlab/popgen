#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)}

library(tools)
library(ggplot2)

pi.fn <- args[1]
my_data <- read.table(pi.fn, sep="\t", header=TRUE)
output_file <- paste(file_path_sans_ext(pi.fn), ".pdf", sep="")

my_data$Chrom <- as.factor(my_data$Chrom)

fig <- ggplot(my_data, aes(x=Chrom, fill=Chrom)) +
  geom_boxplot(aes(ymin=Min., lower=X1st_Qu., middle=Median, upper=X3rd_Qu., ymax=Max.),
               stat = "identity") + ylab("Pi") 

ggsave(output_file, fig, dpi=300, height=5.7, width=15)