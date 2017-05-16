library("pheatmap")
library("data.table")

#load tables into a "list of lists"
qq <- lapply(list.files(".", "*_featurecounts.txt$",full.names=T,recursive=T),function(x) fread(x))
# ensure the samples column is the same name as the treament you want to use:
qq[7]
#merge the "list of lists" into a single table
m <- Reduce(function(...) merge(..., all = T,by=c("Geneid","Chr","Start","End","Strand","Length")), qq)

#convert data.table to data.frame for use with DESeq2
countData <- data.frame(m[,c(1,7:(ncol(m))),with=F])
rownames(countData) <- countData[,1]
countData <- countData[,-1]
indexes <- c("PRO1467_S1","PRO1467_S2","PRO1467_S3", "PRO1467_S10", "PRO1467_S11", "PRO1467_S12", "PRO1467_S4", "PRO1467_S5", "PRO1467_S6", "PRO1467_S14", "PRO1467_S15", "PRO1467_S22", "PRO1467_S7", "PRO1467_S20", "PRO1467_S9", "PRO1467_S17", "PRO1467_S18", "PRO1467_S28")
countDataSubset <- sapply(indexes, function(xx) rowSums(countData[,grep(paste(xx,'_', sep = ""), names(countData)), drop=FALSE]))

#output countData
write.table(countData,"countData.txt",sep="\t",na="",quote=F)
write.table(countDataSubset,"countDataCombined.txt",sep="\t",na="",quote=F)
#output gene details
write.table(m[,1:6,with=F],"genes.txt",sep="\t",quote=F,row.names=F)

library("DESeq2")
library("RColorBrewer")
library("gplots")
library("ggplot2")
library("ggrepel")
#Technical replicates uncollapsed
#orderedColData <- read.table("P.cactorum_RNAseq_design_parsed.txt",header=T,sep="\t")
#colData <- data.frame(orderedColData[ order(orderedColData$Sample.name),])
#unorderedData <- read.table("countData.txt",header=T,sep="\t")
#countData <- data.frame(unorderedData[ , order(colnames(unorderedData))])
#colData$Group <- paste0(colData$Isolate,colData$Plant.Line,colData$Timepoint)

#Technical replicates collapsed
unorderedColData <- read.table("P.cactorum_RNAseq_design_parsed2.txt",header=T,sep="\t")
colData <- data.frame(unorderedColData[ order(unorderedColData$Sample.name),])
unorderedData <- read.table("countDataCombined.txt",header=T,sep="\t")
countData <- data.frame(unorderedData[ , order(colnames(unorderedData))])
colData$Group <- paste0(colData$Isolate,colData$Plant.Line,colData$Timepoint)

design <- ~Group
dds <- DESeqDataSetFromMatrix(countData,colData,design)
sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds, type = c("ratio")))
dds <- DESeq(dds, fitType="local")

vst<-varianceStabilizingTransformation(dds)
pdf("heatmap_vst.pdf", width=12,height=12)
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
pdf("heatmap_rld.pdf")
sampleDists <- dist( t( assay(rld) ) )
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Group)
colnames(sampleDistMatrix) <- paste(rld$Group)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()

###PCA plots
pdf("PCA_vst.pdf")
plotPCA(vst,intgroup=c("Isolate", "Plant.Line", "Timepoint"))
dev.off()

#Plot using rlog transformation:
pdf("PCA_rld.pdf")
plotPCA(rld,intgroup=c("Isolate", "Plant.Line", "Timepoint"))
dev.off()

#Plot using rlog transformation, showing sample names:
data <- plotPCA(rld, intgroup="Group", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

pca_plot<- ggplot(data, aes(PC1, PC2, color=Group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(rld)))
coord_fixed()

ggsave("PCA_sample_names.pdf", pca_plot, dpi=300, height=10, width=12)

###Analysis of gene expression
##"NoEmily0 hours","P414Emily48 hours"
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","NoEmily0 hours","P414Emily48 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"NoEmily0_vs_P414Emily48.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"NoEmily0h_vs_P414Emily48_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"NoEmily0h_vs_P414Emily48_down.txt",sep="\t",na="",quote=F)

##"P414Emily12 hours","P414Emily48 hours"
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","P414Emily12 hours","P414Emily48 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"P414Emily12_vs_P414Emily48.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"P414Emily12_vs_P414Emily48_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"P414Emily12h_vs_P414Emily48_down.txt",sep="\t",na="",quote=F)

##"NoEmily0 hours","P414Emily48 hours"
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","NoEmily0 hours","P414Emily48 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"NoEmily0_vs_P414Emily48.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"NoEmily0_vs_P414Emily48_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"NoEmily0_vs_P414Emily48_down.txt",sep="\t",na="",quote=F)


##"NoFenella0 hours,"P414Fenella48 hours"
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","NoFenella0 hours","P414Fenella48 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"NoFenella0_vs_P414Fenella48.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"NoFenella0h_vs_P414Fenella48_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"NoFenella0h_vs_P414Fenella48_down.txt",sep="\t",na="",quote=F)

##"P414Fenella12 hours","P414Fenella48 hours"
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","P414Fenella12 hours","P414Fenella48 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"P414Fenella12_vs_P414Fenella48.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"P414Fenella12_vs_P414Fenella48_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"P414Fenella12h_vs_P414Fenella48_down.txt",sep="\t",na="",quote=F)

##"NoFenella0 hours","P414Fenella48 hours"
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","NoFenella0 hours","P414Fenella48 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"NoFenella0_vs_P414Fenella48.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"NoFenella0_vs_P414Fenella48_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"NoFenella0_vs_P414Fenella48_down.txt",sep="\t",na="",quote=F)

#"NoFenella0 hours","NoEmily0 hours"
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","NoFenella0 hours","NoEmily0 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"NoFenella0h_vs_NoEmily0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"NoFenella0h_vs_NoEmily0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"NoFenella0h_vs_NoEmily0_down.txt",sep="\t",na="",quote=F)

"P414Fenella12 hours","P414Emily12 hours"
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","P414Fenella12 hours","P414Emily12 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"P414Fenella12h_vs_P414Emily12.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"P414Fenella12h_vs_P414Emily12_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"P414Fenella12h_vs_P414Emily12_down.txt",sep="\t",na="",quote=F)

#"P414Fenella48 hours","P414Emily48 hours"
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","P414Fenella48 hours","P414Emily48 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"P414Fenella48h_vs_P414Emily48.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"P414Fenella48h_vs_P414Emily48_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"P414Fenella48h_vs_P414Emily48_down.txt",sep="\t",na="",quote=F)

#Make a table of raw counts, normalised counts and fpkm values:
raw_counts <- data.frame(counts(dds, normalized=FALSE))
colnames(raw_counts) <- paste(colData$Group)
write.table(raw_counts,"Fvesca_raw_counts2.txt",sep="\t",na="",quote=F)
norm_counts <- data.frame(counts(dds, normalized=TRUE))
colnames(norm_counts) <- paste(colData$Group)
write.table(norm_counts,"Fvesca_normalised_counts.txt",sep="\t",na="",quote=F)

library(Biostrings)
library(naturalsort)
mygenes <- readDNAStringSet("Fragaria_vesca_v1.1.a2_cds.fasta")
t1 <- counts(dds)
t1 <- mygenes[rownames(t1)]
rowRanges(dds) <- GRanges(t1@ranges@NAMES,t1@ranges)

# robust may be better set at false to normalise based on total counts rather than 'library normalisation factors'
fpkm_counts <- data.frame(fpkm(dds, robust = TRUE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"fvesca_fpkm_norm_counts.txt",sep="\t",na="",quote=F)
fpkm_counts <- data.frame(fpkm(dds, robust = FALSE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"fvesca_fpkm_counts.txt",sep="\t",na="",quote=F)