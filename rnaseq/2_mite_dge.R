################################################################################
# 
# Vladimir Zhurov, vzhurov2@uwo.ca, Fri Oct 20 21:43:31 2017
# 
# This is workflow for analysis of RNA-Seq data
#
# Experiment: GENOMITE, mite on tomato data
#
# Interactive analysis using voom/limma
#
# Input:    htseq-count results from tophat/bowtie2 or star
#
# Output:   expression data tables
# 
################################################################################

# Start -------------------------------------------------------------------

# load libraries
library("edgeR")
library("limma")

# specify data set to process
dataset <- "genomite_mite"

# specify mapper used: "star" or "tophat"
mapper <- "star"

# parameters of exactTest and decideTestsDGE
adjust.method = "BH" # default is BH or fdr
p = 0.05
lfc = 1

# set main working path
mainDir <- "/home/sobczm/popgen/rnaseq/genomite/htseq_out/mite_htseq_merged/merged"
dataDir <- "/home/sobczm/popgen/rnaseq/genomite/htseq_out/mite_htseq_merged/merged"
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

# load annotation ---------------------------------------------------------

setwd(annoDir)
annotation.length <- read.delim("tetur_transcript_length-current.txt", header = TRUE)
annotation <- read.delim("tetur_annotation-current.txt", header = TRUE)


# voom analysis of differential gene expression --------------------------

# load count data
setwd(dataDir)

# load Targets.txt file with file list, groups and labels
targets <- readTargets()
dge <- readDGE(files = targets$File, header = FALSE)

# calculate library size normalization factors
dge <- calcNormFactors(dge, method = "TMM")

# keep loci that are expressed
isexpr <- rowSums(cpm(dge) > 1) >= 3
dge <- dge[isexpr, keep.lib.sizes = FALSE]

# Constitutive differences between mite populations -----------------------

# create design matrix
batch <- factor(targets$Batch)
treatment <- factor(targets$Mite, levels = c("nonadapted", "adapted"))
design <- model.matrix(~ batch + treatment)
colnames(design) <- c("Intercept", "batch.2", "batch.3", "adapt")

#  apply voom transformation
# dge.voom <- voom(dge, design, plot = FALSE)
dge.voom <- voomWithQualityWeights(dge, design, plot = FALSE)

# fit model
fit <- lmFit(dge.voom, design)

# set suffix for output files
lfc_text <- round(2^lfc, digits = 1)
suffix <- paste("_voomW_Mite", mapper, lfc_text, p, adjust.method, sep = "_")

# create and fit contrasts
cont <- makeContrasts(
  Mite.Adapt.v.Nonadapt = adapt,
  levels = design
)
c.fit <- contrasts.fit(fit, cont)
eb.fit <- eBayes(c.fit)

### generate topTable based on any time point and filter using decideTests

topTable.complete <- topTable(eb.fit, number = Inf, sort.by = "none", adjust.method = adjust.method)
topTable.complete <- cbind(ID = rownames(topTable.complete), topTable.complete)
topTable.complete <- merge(topTable.complete, annotation, by = "ID", all.x = TRUE)
# apply cut-off's
dt <- decideTests(eb.fit, lfc = lfc, p = p, adjust.method = adjust.method, method = "separate")
rs <- rowSums(abs(dt))
topTable.filtered <- topTable.complete[rs != 0, ]
# save limma result
setwd(limmaDir)
filename <- paste("results_", dataset, suffix, "_all-individual_p_based.txt", sep = "")
write.table(x = topTable.filtered, file = filename, row.names = FALSE, sep = "\t", quote = FALSE)
# save gene list
setwd(geneListsDir)
filename.list <- paste("list_", dataset, suffix, "_all-individual_p_based.txt", sep = "")
write.table(x = topTable.filtered$ID, file = filename.list, quote = FALSE, col.names = FALSE, row.names = FALSE)

### generate DEG table based on overall p-value

topTable.overall <- topTable(eb.fit, number = Inf, p = p, lfc = lfc, adjust.method = adjust.method)
topTable.overall <- cbind(ID = rownames(topTable.overall), topTable.overall)
topTable.overall <- merge(topTable.overall, annotation, by = "ID", all.x = TRUE)
# save limma result
setwd(limmaDir)
filename <- paste("results_", dataset, suffix, "_all-overall_p_based.txt", sep = "")
write.table(x = topTable.overall, file = filename, row.names = FALSE, sep = "\t", quote = FALSE)
# save gene list
setwd(geneListsDir)
filename.list <- paste("list_", dataset, suffix, "_all-overall_p_based.txt", sep = "")
write.table(x = topTable.overall$ID, file = filename.list, quote = FALSE, col.names = FALSE, row.names = FALSE)

### save individual results and gene lists per contrast

l <- dim(cont)
for (i in 1:l[2]) {
  # generate topTable for a given contrast
  topTable <- topTable(eb.fit, coef = i, genelist = fit$genes, n = Inf, p = p, lfc = lfc, adjust.method = adjust.method)
  topTable <- cbind(ID = rownames(topTable), topTable)
  topTable <- merge(topTable, annotation, by = "ID", all.x = TRUE)
  cont.names <- colnames(cont)
  # save limma result
  setwd(limmaDir)
  filename.result <- paste0("results_", dataset, suffix, "_", cont.names[i], ".txt")
  write.table(x = topTable, file = filename.result, row.names = FALSE, sep = "\t", quote = FALSE)
  # save gene list
  setwd(geneListsDir)
  filename.list <- paste0("list_", dataset, suffix, "_", cont.names[i], ".txt")
  write.table(x = topTable$ID, file = filename.list, quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  # separate up- and down-regulated genes
  
  # upregulated genes
  topTable.up <- topTable[(topTable$logFC > 0), ]
  # save results
  setwd(limmaDir)
  filename.result <- paste0("results_", dataset, suffix, "_", cont.names[i], "-Up.txt")
  write.table(x = topTable.up, file = filename.result, row.names = FALSE, sep = "\t")
  # save gene list
  setwd(geneListsDir)
  filename.list <- paste0("list_", dataset, suffix, "_", cont.names[i], "-Up.txt")
  write.table(x = topTable.up$ID, file = filename.list, quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  # downregulated genes
  topTable.down <- topTable[(topTable$logFC < 0), ]
  # save results
  setwd(limmaDir)
  filename.result <- paste0("results_", dataset, suffix, "_", cont.names[i], "-Down.txt")
  write.table(x = topTable.down, file = filename.result, row.names = FALSE, sep = "\t")
  # save gene list
  setwd(geneListsDir)
  filename.list <- paste0("list_", dataset, suffix, "_", cont.names[i], "-Down.txt")
  write.table(x = topTable.down$ID, file = filename.list, quote = FALSE, col.names = FALSE, row.names = FALSE) 
}


# Constitutive difference between drought and no stress --------------------

# create design matrix
batch <- factor(targets$Batch)
treatment <- factor(targets$Plant, levels = c("cntrl", "nostress", "drought"))
design <- model.matrix(~ batch + treatment)
colnames(design) <- c("Intercept", "batch.2", "batch.3", "nodrought", "drought")

#  apply voom transformation
# dge.voom <- voom(dge, design, plot = FALSE)
dge.voom <- voomWithQualityWeights(dge, design, plot = FALSE)

# fit model
fit <- lmFit(dge.voom, design)

# set LFC and FDR cut-offs
p = 0.05
lfc = 0.58 # absolute fold change of 1.5, 50% increase or decrease

# set suffix for output files
lfc_text <- round(2^lfc, digits = 1)
suffix <- paste("_voomW_Drought", mapper, lfc_text, p, adjust.method, sep = "_")

# create and fit contrasts
cont <- makeContrasts(
  Nodrought.Cntrl = nodrought,
  Drought.Cntrl = drought,
  Drought.Nodrought = drought - nodrought,
  levels = design
)
c.fit <- contrasts.fit(fit, cont)
eb.fit <- eBayes(c.fit)

### generate topTable based on any time point and filter using decideTests

topTable.complete <- topTable(eb.fit, number = Inf, sort.by = "none", adjust.method = adjust.method)
topTable.complete <- cbind(ID = rownames(topTable.complete), topTable.complete)
topTable.complete <- merge(topTable.complete, annotation, by = "ID", all.x = TRUE)
# apply cut-off's
dt <- decideTests(eb.fit, lfc = lfc, p = p, adjust.method = adjust.method, method = "separate")
rs <- rowSums(abs(dt))
topTable.filtered <- topTable.complete[rs != 0, ]
# save limma result
setwd(limmaDir)
filename <- paste("results_", dataset, suffix, "_all-individual_p_based.txt", sep = "")
write.table(x = topTable.filtered, file = filename, row.names = FALSE, sep = "\t", quote = FALSE)
# save gene list
setwd(geneListsDir)
filename.list <- paste("list_", dataset, suffix, "_all-individual_p_based.txt", sep = "")
write.table(x = topTable.filtered$ID, file = filename.list, quote = FALSE, col.names = FALSE, row.names = FALSE)

### generate DEG table based on overall p-value

topTable.overall <- topTable(eb.fit, number = Inf, p = p, lfc = lfc, adjust.method = adjust.method)
topTable.overall <- cbind(ID = rownames(topTable.overall), topTable.overall)
topTable.overall <- merge(topTable.overall, annotation, by = "ID", all.x = TRUE)
# save limma result
setwd(limmaDir)
filename <- paste("results_", dataset, suffix, "_all-overall_p_based.txt", sep = "")
write.table(x = topTable.overall, file = filename, row.names = FALSE, sep = "\t", quote = FALSE)
# save gene list
setwd(geneListsDir)
filename.list <- paste("list_", dataset, suffix, "_all-overall_p_based.txt", sep = "")
write.table(x = topTable.overall$ID, file = filename.list, quote = FALSE, col.names = FALSE, row.names = FALSE)

### save individual results and gene lists per contrast

l <- dim(cont)
for (i in 1:l[2]) {
  # generate topTable for a given contrast
  topTable <- topTable(eb.fit, coef = i, genelist = fit$genes, n = Inf, p = p, lfc = lfc, adjust.method = adjust.method)
  topTable <- cbind(ID = rownames(topTable), topTable)
  topTable <- merge(topTable, annotation, by = "ID", all.x = TRUE)
  cont.names <- colnames(cont)
  # save limma result
  setwd(limmaDir)
  filename.result <- paste0("results_", dataset, suffix, "_", cont.names[i], ".txt")
  write.table(x = topTable, file = filename.result, row.names = FALSE, sep = "\t", quote = FALSE)
  # save gene list
  setwd(geneListsDir)
  filename.list <- paste0("list_", dataset, suffix, "_", cont.names[i], ".txt")
  write.table(x = topTable$ID, file = filename.list, quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  # separate up- and down-regulated genes
  
  # upregulated genes
  topTable.up <- topTable[(topTable$logFC > 0), ]
  # save results
  setwd(limmaDir)
  filename.result <- paste0("results_", dataset, suffix, "_", cont.names[i], "-Up.txt")
  write.table(x = topTable.up, file = filename.result, row.names = FALSE, sep = "\t")
  # save gene list
  setwd(geneListsDir)
  filename.list <- paste0("list_", dataset, suffix, "_", cont.names[i], "-Up.txt")
  write.table(x = topTable.up$ID, file = filename.list, quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  # downregulated genes
  topTable.down <- topTable[(topTable$logFC < 0), ]
  # save results
  setwd(limmaDir)
  filename.result <- paste0("results_", dataset, suffix, "_", cont.names[i], "-Down.txt")
  write.table(x = topTable.down, file = filename.result, row.names = FALSE, sep = "\t")
  # save gene list
  setwd(geneListsDir)
  filename.list <- paste0("list_", dataset, suffix, "_", cont.names[i], "-Down.txt")
  write.table(x = topTable.down$ID, file = filename.list, quote = FALSE, col.names = FALSE, row.names = FALSE) 
}

# End ---------------------------------------------------------------------

