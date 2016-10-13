#!/usr/bin/env Rscript

#Print a bootstrapped (100x) Newick tree based on input FASTA alignment.
#First argument (required): input FASTA alignment file

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)}


library(ape)
library(tools)

input <- args[1]
out <- read.dna(input, format="fasta")
f <- function(x) njs(dist.dna(x))
#Construct initial tree
tree <- f(out)
#100 bootrstrap replicates
bstrees <- boot.phylo(tree, out, f, trees = TRUE)$trees
## get proportions of each bipartition and draw the tree, printing it to a PDF file
boot <- prop.clades(tree, bstrees)
tree_pdf <- paste(file_path_sans_ext(input), "_nj.pdf", sep="")
pdf(tree_pdf, width=6, height=8)
plot(tree)
drawSupportOnEdges(boot)
dev.off()
#Save the bootstrap values in the PHYLO object
tree$node.label <- boot
#Write Newick tree
tree_out <- paste(file_path_sans_ext(input), "_nj.nwk", sep="")
write.tree(tree, file=tree_out)
