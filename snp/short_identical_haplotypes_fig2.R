#!/usr/bin/env Rscript
library("ggplot2")
library("tools")
library("plyr")
#Plot mean percentage perfect matching fragments per sample

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
plot.data.diff2 <- subset(stats, perc_matching_haplo != "NA")
res.plyr <- ddply(plot.data.diff2, .(sample_name), function(x) mean(x$perc_matching_haplo) )
fig <- ggplot(res.plyr, aes(x=sample_name, y=1, fill=sample_name)) + geom_bar(position=position_dodge(), stat="identity") +
  ylab("% matching haplotype fragments") + scale_y_continuous(expand=c(0,0)) + theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())
output_file <- paste(file_path_sans_ext(input_file), "_2.pdf", sep="")
ggsave(output_file, fig, dpi=300, width=as.numeric(args[2]), height=as.numeric(args[3]))

