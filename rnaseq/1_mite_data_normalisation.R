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

# keep unfitered by expressison DGE object
dge.all <- dge

# keep loci that are expressed
isexpr <- rowSums(cpm(dge) > 1) >= 3
dge <- dge[isexpr, keep.lib.sizes = FALSE]

# Save normalized count data ----------------------------------------------

# CPM data

# save normalized CPM counts
norm.count <- cpm(dge)
norm.count.table <- cbind(ID = rownames(norm.count), norm.count)
setwd(cpmDir)
norm.count.file <- paste0("normalized_gene_expression-CPM-", dataset, ".txt")
write.table(norm.count.table, file = norm.count.file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# log-transformed CPM counts
norm.count.log <- cpm(dge, log = TRUE)
norm.count.log.table <- cbind(ID = rownames(norm.count.log), norm.count.log)
setwd(cpmDir)
norm.count.file <- paste0("normalized_gene_expression-log2CPM-", dataset, ".txt")
write.table(norm.count.log.table, file = norm.count.file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# remove batch effect from log-transformed CPM counts
lev <- unique(targets$Batch)
batch <- factor(targets$Batch, levels = lev)
norm.count.log <- removeBatchEffect(norm.count.log, batch = batch)
norm.count.log.table <- cbind(ID = rownames(norm.count.log), norm.count.log)
setwd(cpmDir)
norm.count.file <- paste0("normalized_gene_expression-removed_batch_effect-log2CPM-", dataset, ".txt")
write.table(norm.count.log.table, file = norm.count.file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# FPKM data

# get gene length
count.table <- dge$counts
count.table <- cbind("ID" = row.names(count.table), count.table)
count.table <- merge(count.table, annotation.length[, c(1, 2)], by = "ID", all.x = TRUE, sort = FALSE)
gene.length <- as.vector(count.table$Length)
gene.length <- as.numeric(gene.length)
# calculate FPKM
rpkm.count <- rpkm(dge, gene.length = gene.length)
rpkm.count.table <- cbind(ID = rownames(rpkm.count), rpkm.count)
setwd(cpmDir)
rpkm.count.file <- paste0("normalized_gene_expression-FPKM-", dataset, ".txt")
write.table(rpkm.count.table, file = rpkm.count.file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# log-transformed FPKM counts
rpkm.count.log <- rpkm(dge, gene.length = gene.length, log = TRUE)
rpkm.count.log.table <- cbind(ID = rownames(rpkm.count.log), rpkm.count.log)
setwd(cpmDir)
rpkm.count.file <- paste0("normalized_gene_expression-log2FPKM-", dataset, ".txt")
write.table(rpkm.count.log.table, file = rpkm.count.file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# remove batch effect from log-transformed FPKM counts
lev <- unique(targets$Batch)
batch <- factor(targets$Batch, levels = lev)
rpkm.count.log <- removeBatchEffect(rpkm.count.log, batch = batch)
rpkm.count.log.table <- cbind(ID = rownames(rpkm.count.log), rpkm.count.log)
setwd(cpmDir)
rpkm.count.file <- paste0("normalized_gene_expression-removed_batch_effect-log2FPKM-", dataset, ".txt")
write.table(rpkm.count.log.table, file = rpkm.count.file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# Save unfiltered normalized count data -----------------------------------

# CPM data

# save normalized CPM counts
norm.count <- cpm(dge.all)
norm.count.table <- cbind(ID = rownames(norm.count), norm.count)
setwd(cpmDir)
norm.count.file <- paste0("normalized_all_gene_expression-CPM-", dataset, ".txt")
write.table(norm.count.table, file = norm.count.file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# log-transformed CPM counts
norm.count.log <- cpm(dge.all, log = TRUE)
norm.count.log.table <- cbind(ID = rownames(norm.count.log), norm.count.log)
setwd(cpmDir)
norm.count.file <- paste0("normalized_all_gene_expression-log2CPM-", dataset, ".txt")
write.table(norm.count.log.table, file = norm.count.file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# remove batch effect from log-transformed CPM counts
lev <- unique(targets$Batch)
batch <- factor(targets$Batch, levels = lev)
norm.count.log <- removeBatchEffect(norm.count.log, batch = batch)
norm.count.log.table <- cbind(ID = rownames(norm.count.log), norm.count.log)
setwd(cpmDir)
norm.count.file <- paste0("normalized_all_gene_expression-removed_batch_effect-log2CPM-", dataset, ".txt")
write.table(norm.count.log.table, file = norm.count.file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# FPKM data

# get gene length
count.table <- dge.all$counts
count.table <- cbind("ID" = row.names(count.table), count.table)
count.table <- merge(count.table, annotation.length[, c(1, 2)], by = "ID", all.x = TRUE, sort = FALSE)
gene.length <- as.vector(count.table$Length)
gene.length <- as.numeric(gene.length)
# calculate FPKM
rpkm.count <- rpkm(dge.all, gene.length = gene.length)
rpkm.count.table <- cbind(ID = rownames(rpkm.count), rpkm.count)
setwd(cpmDir)
rpkm.count.file <- paste0("normalized_all_gene_expression-FPKM-", dataset, ".txt")
write.table(rpkm.count.table, file = rpkm.count.file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

rpkm.count.log <- rpkm(dge.all, gene.length = gene.length, log = TRUE)
rpkm.count.log.table <- cbind(ID = rownames(rpkm.count.log), rpkm.count.log)
setwd(cpmDir)
rpkm.count.file <- paste0("normalized_all_gene_expression-log2FPKM-", dataset, ".txt")
write.table(rpkm.count.log.table, file = rpkm.count.file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# remove batch effect from log-transformed FPKM counts
lev <- unique(targets$Batch)
batch <- factor(targets$Batch, levels = lev)
rpkm.count.log <- removeBatchEffect(rpkm.count.log, batch = batch)
rpkm.count.log.table <- cbind(ID = rownames(rpkm.count.log), rpkm.count.log)
setwd(cpmDir)
rpkm.count.file <- paste0("normalized_all_gene_expression-removed_batch_effect-log2FPKM-", dataset, ".txt")
write.table(rpkm.count.log.table, file = rpkm.count.file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# End ---------------------------------------------------------------------
