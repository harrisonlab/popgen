#!/usr/bin/env Rscript
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
plot.data <- data.frame(sample_name=stats$sample_name, perc_matching_haplo=stats$perc_matching_haplo)
fig <- ggplot(plot.data, aes(x=sample_name, y=perc_matching_haplo, fill=sample_name)) + geom_bar(position=position_dodge(), stat="identity") + ylab("% haplotype match") + scale_y_continuous(expand=c(0,0)) + theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())
output_file <- paste(file_path_sans_ext(input_file), ".pdf", sep="")
ggsave(output_file, fig, dpi=300, width=as.numeric(args[2]), height=as.numeric(args[3]))

