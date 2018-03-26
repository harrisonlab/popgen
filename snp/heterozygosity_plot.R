#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)}
library(tools)
library(ggplot2)

het.fn <- args[1]
my_data <- read.table(het.fn, sep="\t", header=TRUE)
output_file <- paste(file_path_sans_ext(het.fn), ".pdf", sep="")

my_data <- my_data[order(my_data$F),]
# lock in factor level order
my_data$INDV <- factor(my_data$INDV, levels=my_data$INDV)

fig <- 

ggplot(my_data, aes(x=INDV, y=my_data$F)) + geom_bar(position=position_dodge(), stat="identity") +
  coord_flip() + ylab("F het")

ggsave(output_file, fig, dpi=300, height=40, width=6.13)