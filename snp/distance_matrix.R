#!/usr/bin/env Rscript

#Print a heatmap displaynig the distance results output by similarity_percentage.py
#First argument (required): input distance results
#Second argument (optional): output filename
library(tools)
args = commandArgs(trailingOnly=TRUE)


if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = paste(file_path_sans_ext(args[1]), "_distance_matrix.pdf", sep="")
}

library(gplots)
library(RColorBrewer)
distance = read.table(args[1], header=T, stringsAsFactors=F)
#Create basic palette
mypalette<-brewer.pal(9, "YlGn")
#Ramp up the color palette
cols <- colorRampPalette (mypalette) (100)
par(mar=c(9,6,6,4)+0.1) 
pdf(args[2], width=200, height=130)
heatmap.2(as.matrix(distance),dendrogram="col", cellnote=as.matrix(distance),
          notecol="black",col=cols,scale="none",key=TRUE, keysize=0.25,
          density.info="none", cexRow=2.0,cexCol=2.0, notecex=2.0, trace="none",margins=c(30,30))
dev.off()
