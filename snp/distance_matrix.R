#!/usr/bin/env Rscript

#Print a heatmap displaynig the distance results output by similarity_percentage.py
#First argument (required): input distance results
#Second argument (optional): output PDF filename

args = commandArgs(trailingOnly=TRUE)


if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "distance_matrix.pdf"
}

library(gplots)
library(RColorBrewer)
distance = read.table(args[1], header=T, stringsAsFactors=F)
#Create basic palette
mypalette<-brewer.pal(9, "YlGn")
#Ramp up the color palette
cols <- colorRampPalette (mypalette) (100)
pdf(args[2], width=11, height=8)
heatmap.2(as.matrix(distance),dendrogram="col", cellnote=as.matrix(distance),
          notecol="black",col=cols,scale="none",key=TRUE, keysize=1.5,
          density.info="none", trace="none", cexRow=1.0,cexCol=1.0)
dev.off()
