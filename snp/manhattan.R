library(qqman)
library(tools)
args = commandArgs(trailingOnly=TRUE)
input = args[1]
output <- paste(input, "_manhattan.pdf", sep="")
stats <- read.csv(input, sep="\t", quote='', stringsAsFactors=TRUE,header=TRUE)
pdf(output,7,5)
manhattan(stats, col = c("blue", "red3"), ymax = 12)
dev.off()