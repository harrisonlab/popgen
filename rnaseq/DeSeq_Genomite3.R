library("pheatmap")
library("data.table")

#load tables into a "list of lists"
qq <- lapply(list.files(".", "*_fc.out$",full.names=T,recursive=T),function(x) fread(x))
# ensure the samples column is the same name as the treament you want to use:
qq[7]
#merge the "list of lists" into a single table
m <- Reduce(function(...) merge(..., all = T,by=c("Geneid","Chr","Start","End","Strand","Length")), qq)
sample_names <- colnames(m)[7:78]

#convert data.table to data.frame for use with DESeq2
countData <- data.frame(m[,c(1,7:(ncol(m))),with=F])
rownames(countData) <- countData[,1]
countData <- countData[,-1]

write.table(countData,"strawberry_countData.txt",sep="\t",na="",quote=F)
write.table(m[,1:6,with=F],"genes.txt",sep="\t",quote=F,row.names=F)

library("DESeq2")
library("RColorBrewer")
library("gplots")
library("ggplot2")
library("ggrepel")
#Technical replicates uncollapsed
orderedColData <- read.table("strawberry_deseq_samples.txt",header=T,sep="\t")
colData <- data.frame(orderedColData[ order(orderedColData$Sample.name),])
unorderedData <- read.table("strawberry_countData.txt",header=T,sep="\t")
countData <- data.frame(unorderedData[ , order(colnames(unorderedData))])
colData$Group <- paste0(colData$Stress, colData$Mite, colData$Time)
design <- ~Group
dds <- DESeqDataSetFromMatrix(countData,colData,design)
sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds, type = c("ratio")))
dds <- DESeq(dds, fitType="local")

vst<-varianceStabilizingTransformation(dds)
pdf("strawberry_heatmap_vst.pdf", width=36,height=36)
sampleDists<-dist(t(assay(vst)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$Group)
colnames(sampleDistMatrix) <- paste(vst$Group)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix,
         trace="none",  # turns off trace lines inside the heat map
         col=colours, # use on color palette defined earlier
         margins=c(12,12), # widens margins around plot
         srtCol=45,
         srtCol=45)
dev.off()

# Sample distances measured with rlog transformation:
rld <- rlog( dds )
pdf("strawberry_heatmap_rld.pdf", width=36,height=36)
sampleDists <- dist( t( assay(rld) ) )
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Group)
colnames(sampleDistMatrix) <- paste(rld$Group)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(36,36),srtCol=45)
dev.off()

##PCA plots
pdf("strawberry_PCA_vst.pdf", width=12,height=12)
plotPCA(vst,intgroup=c("Stress", "Mite", "Time"))
dev.off()

#Plot using rlog transformation:
pdf("strawberry_PCA_rld.pdf", width=12,height=12)
plotPCA(rld,intgroup=c("Stress", "Mite", "Time"))
dev.off()

#Plot using rlog transformation, showing sample names:
data <- plotPCA(rld, intgroup="Group", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

pca_plot<- ggplot(data, aes(PC1, PC2, color=Group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(rld)))
coord_fixed()

ggsave("strawberry_PCA_sample_names.pdf", pca_plot, dpi=300, height=12, width=12)
