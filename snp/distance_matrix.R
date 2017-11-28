#!/usr/bin/env Rscript

#Print a heatmap displaynig the distance results output by similarity_percentage.py
#First argument (required): input distance results
#Second argument (optional) - output figure width
#Third argument (optional) - output figure height

###Note for usage:
##This script takes a bit of work to produce the plot matching the size of your dataset. 
##The following parameters inside the script may have to be changed:
#Margins - Specify the margins size before running the heatmap.2 command
#cexRow - Row label size
#cexCol - Column label size
#notecex - Size of values within individual cells
#INTeger in colorRampPalette (mypalette) (INT) - depending on the range of the values of % match, may want to increase or decrease number of colors in the scale as specified here.

library(tools)
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = 16
  args[3] = 8
}

output <- paste(file_path_sans_ext(args[1]), "_distance_matrix.pdf", sep="")

library(gplots)
library(RColorBrewer)
distance = read.table(args[1], header=T, stringsAsFactors=F)
#Create basic palette
mypalette<-brewer.pal(9, "YlGn")
#Ramp up the color palette
cols <- colorRampPalette (mypalette) (100)
par(mar=c(30,5,5,30)+0.1) 
pdf(output, as.numeric(args[2]), as.numeric(args[3]))
heatmap.2(as.matrix(distance),dendrogram="col", cellnote=as.matrix(distance),
          notecol="black",col=cols,scale="none",key=TRUE, keysize=0.75, margins=c(15,15),
          density.info="none", trace="none", cexRow=1, cexCol=1, notecex=0.7)
dev.off()
