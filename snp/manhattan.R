library(qqman)
library(tools)
args = commandArgs(trailingOnly=TRUE)
input = args[1]
output <- paste(input, "_manhattan.pdf", sep="")
stats <- read.csv(input, sep="\t", quote='', stringsAsFactors=TRUE,header=TRUE)
stats2 <- stats[!is.na(stats$P),]
pdf(output,30,5)
manhattan(stats2, col = c("blue", "red3"), ymax = 12)
dev.off()