#!/usr/bin/env Rscript
#Scatterplot of intervals showing non-perfect matches
library("ggplot2")
library("tools")
#First argument - stats filename 
#Second argument (optional) - output figure width
#Third argument (optional) - output figure height

args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) 
{
  stop("You must supply at least one argument.")
} else if (length(args)==1) 
{
  args[2] = 8
  args[3] = 5.7
}

input_file <- args[1]
stats <- read.csv(input_file, sep="\t", quote='', stringsAsFactors=TRUE,header=TRUE)
#Remove 100% and NA matches.
plot.data.diff <- subset(stats, perc_matching_haplo != 100 & perc_matching_haplo != "NA")
#Scatterplot of intervals showing non-perfect matches
fig <- ggplot(plot.data.diff, aes(x=interval_start,y=perc_matching_haplo,color=sample_name)) + geom_point() + ylab("% identity of haplotype fragments") + scale_x_continuous(breaks = pretty(plot.data.diff$interval_start, n = 20)) 
output_file <- paste(file_path_sans_ext(input_file), "_1.pdf", sep="")
ggsave(output_file, fig, dpi=300, width=as.numeric(args[2]), height=as.numeric(args[3]))
