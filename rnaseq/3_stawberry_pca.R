################################################################################
# 
# Vladimir Zhurov, vzhurov2@uwo.ca, Tue Aug  1 11:03:59 2017
# 
# This is workflow for analysis of RNA-Seq data
#
# Experiment: GENOMITE, mite on tomato data
#
# This is the sub-routine for PCA
#
# Step 01:  generate 2D PCA plot
#
# Step 02:  generate interactive 3D PCA plot and animation
#
# Input:    voom processed dge and topTable results
#
# Output:   2D and 3D PCA plots and animations
#
################################################################################

library("limma")


# settings ----------------------------------------------------------------

# specify data set to process
dataset <- "genomite_strawberry"

# specify mapper used: "star" or "tophat"
mapper <- "star"

# parameters of exactTest and decideTestsDGE
adjust.method = "BH" # default is BH or fdr
p = 0.05
lfc = 1

# set main working path
mainDir <- "/home/sobczm/popgen/rnaseq/genomite/htseq_out/strawberry_htseq_merged/merged"
dataDir <- "/home/sobczm/popgen/rnaseq/genomite/htseq_out/strawberry_htseq_merged/merged"
annoDir <- paste0(mainDir, "/annotation")
edgerDir <- paste0(mainDir, "/edger")
voomDir <- paste0(mainDir, "/voom")
limmaDir <- voomDir
cpmDir <- paste0(mainDir, "/cpm")
geneListsDir <- paste0(mainDir, "/genelists")
gsaDir <- paste0(mainDir, "/gsa")
graphsDir <- paste0(mainDir, "/graphs")
topGODir <- paste0(mainDir, "/topGO")
heatmapDir <- paste0(mainDir, "/heatmap")
pcaDir <- paste0(mainDir, "/pca")

# functions ---------------------------------------------------------------

color.map <- function(groups) {
  color.number <- length(unique(groups))
  colors <- rainbow(color.number)
  group.color <- colors[groups]
  return(group.color)
}

# main analysis -----------------------------------------------------------
setwd(dataDir)
data.original <- rpkm.count.log
data <- data.original

# remove batch effect 
lev <- unique(targets$Batch)
batch <- factor(targets$Batch, levels = lev)
data <- removeBatchEffect(data, batch = batch)

# perform PCA
pca.result <- prcomp(t(data), retx = TRUE, center = TRUE, scale = FALSE)
perc.expl <- round(((pca.result$sdev[1:3]^2)/sum(pca.result$sdev^2))*100, 1)

# what to color samples by
# line color
# targets <- targets[-c(1:6, 13:23, 35:45), ]
color.by <- targets$Time
color.by <- as.numeric(factor(color.by))
line.color <- color.map(color.by)
# fill color
# option 1
fill.color <- "white"
# option 2
color.by <- targets$Time
color.by <- as.numeric(factor(color.by))
fill.color <- color.map(color.by)

# 2D plot -----------------------------------------------------------------

par(mfrow = c(2, 2), oma = c(0, 0, 1, 0))

plot(pca.result$x[, 1], pca.result$x[, 2], pch = 16, col = fill.color,
     xlab = paste("PC1 (", perc.expl[1], "%)",sep = ""),
     ylab = paste("PC2 (", perc.expl[2], "%)",sep = ""))
text(pca.result$x[, 1], pca.result$x[, 2], targets$Label, pos = 1, cex = 0.75)
points(pca.result$x[, 1], pca.result$x[, 2], pch = 1, col = line.color)

plot(pca.result$x[, 1], pca.result$x[, 3], pch = 16, col = fill.color,
     xlab = paste("PC1 (", perc.expl[1], "%)",sep = ""),
     ylab = paste("PC3 (", perc.expl[3], "%)",sep = ""))
text(pca.result$x[, 1], pca.result$x[, 3], targets$Label, pos = 1, cex = 0.75)
points(pca.result$x[, 1], pca.result$x[, 3], pch = 1, col = line.color)

plot(pca.result$x[, 2], pca.result$x[, 3], pch = 16, col = fill.color,
     xlab = paste("PC2 (", perc.expl[2], "%)",sep = ""),
     ylab = paste("PC3 (", perc.expl[3], "%)",sep = ""))
text(pca.result$x[, 2], pca.result$x[, 3], targets$Label, pos = 1, cex = 0.75)
points(pca.result$x[, 2], pca.result$x[, 3], pch = 1, col = line.color)

explained.by.three <- perc.expl[1] + perc.expl[2] + perc.expl[3]

barplot((100*pca.result$sdev^2)/sum(pca.result$sdev^2), ylim = c(0, 80),
        xlab = paste("Principal components\nFirst three PCs: ", explained.by.three, "%", sep =""),
        ylab = "% of total variance explained")

title("Color: Time status", outer = TRUE)

par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))


