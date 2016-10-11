library("PopGenome")
library(ggplot2)

######################## BEFORE RUNNING ################################
#Assign individuals to appropriate populations (or just 1!)
#This script calculates: Fst all/pairwise per gene or sliding window or per SNP
#More than one population needs to be defined, of course!
nonpatho <- c("FOCA1-2", "FOCA28", "FOCCB3", "FOCD2", "FOCHB6", "FOCPG")
patho <- c("FOCA23", "FOC55", "FOC125", "FOCFus2")
#Need to set argument diploid=TRUE if using diploid genomes in the below command:
populations <- list(nonpatho, patho)
#Number of populations assigned above.
population_no <- length(populations)
population_names <- c("nonpatho", "patho") #Given in the same order, as above.
#Interval and jump size used in the sliding window analysis
interval <-  1000
jump_size <-  interval / 10
#########################################################################

#Folder containing FASTA alignments in the current dir
gff <- "gff"
all_folders <- list.dirs("contigs", full.names = FALSE)
#Remove the gff folder from PopGenome contig analysis
contig_folders <- all_folders[all_folders != "gff"]

###Loop through each contig-containing folder to calculate stats on each contig separately.
for (dir in contig_folders[contig_folders != ""])
{
  contig_folder <- paste("contigs/", dir, sep="")
  GENOME.class <- readData(contig_folder, gffpath=gff, include.unknown = TRUE)
  GENOME.class <- set.populations(GENOME.class, populations)

###############################################################################

#### Gene-based analysis
GENOME.class.split <- splitting.data(GENOME.class, subsites="gene")
GENOME.class.split <- F_ST.stats(GENOME.class.split)
get.F_ST(GENOME.class.split)

FST_all <- GENOME.class.split@nuc.F_ST.vs.all
FST_all_d <- as.data.frame(FST_all)
FST_pairwise <- GENOME.class.split@nuc.F_ST.pairwise
FST_pairwise_d <- as.data.frame(as.vector(FST_pairwise))
Hudson_KST <- GENOME.class.split@Hudson.K_ST
Hudson_KST_d <- as.data.frame(Hudson_KST)

for (i in seq_along(population_names))
{
  file_hist <- paste(dir, "_", population_names[i], "_total_FST_per_gene.pdf", sep="")
  fst_plot <- ggplot(FST_all_d, aes(x=FST_all_d[,i])) + geom_histogram(colour="black", fill="darkseagreen") + ggtitle(dir) + xlab(expression(paste("Total FST per gene"))) + ylab("Number of genes") + scale_x_continuous(breaks = pretty(FST_all_d[,i], n = 10))
  ggsave(file_hist, fst_plot)
  file_table = paste(dir, "_", population_names[i], "_total_FST_per_gene.txt", sep="")
  current_gff <- paste(gff, "/", dir, ".gff", sep="")
  gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr=dir, feature=FALSE, extract.gene.names=TRUE)
  fst_table <- cbind(gene_ids, FST_all[,i])
  write.table(fst_table, file=file_table, sep="\t",quote=FALSE, col.names=FALSE)
}

for (i in seq_along(nrow(FST_pairwise)))
{
  file_hist <- paste(dir, "_", population_names[i], "_pairwise_FST_per_gene", i, ".pdf", sep="")
  fst_plot <- ggplot(FST_pairwise_d, aes(x=FST_pairwise_d[,i])) + geom_histogram(colour="black", fill="cadetblue") + ggtitle(dir) + xlab(expression(paste("Pairwise FST per gene"))) + ylab("Number of genes") + scale_x_continuous(breaks = pretty(FST_pairwise_d[,i], n = 10))
  ggsave(file_hist, fst_plot)
  file_table = paste(dir, "_", population_names[i], "_pairwise_FST_per_gene", i, ".txt", sep="")
  current_gff <- paste(gff, "/", dir, ".gff", sep="")
  gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr=dir, feature=FALSE, extract.gene.names=TRUE)
  fst_table <- cbind(gene_ids, FST_pairwise_d[,i])
  write.table(fst_table, file=file_table, sep="\t",quote=FALSE, col.names=FALSE, row.names=FALSE)
}

### Plot Hudson KST

for (i in seq_along(Hudson_KST_d))
{

  file_hist <- paste(dir, "_", population_names[i], "_Hudson_KST_per_gene", i, ".pdf", sep="")
  fst_plot <- ggplot(Hudson_KST_d, aes(x=Hudson_KST_d[,i])) + geom_histogram(colour="black", fill="springgreen") + ggtitle(dir) + xlab(expression(paste("Hudson KST per gene"))) + ylab("Number of genes") + scale_x_continuous(breaks = pretty(Hudson_KST_d[,i], n = 10))
  ggsave(file_hist, fst_plot)
  file_table = paste(dir, "_", population_names[i], "_Hudson_KST_per_gene", i, ".txt", sep="")
  current_gff <- paste(gff, "/", dir, ".gff", sep="")
  gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr=dir, feature=FALSE, extract.gene.names=TRUE)
  fst_table <- cbind(gene_ids, Hudson_KST_d[,i])
  write.table(fst_table, file=file_table, sep="\t",quote=FALSE, col.names=FALSE, row.names=FALSE)

}

#### Sliding window analysis (interval)
GENOME.class.slide <- sliding.window.transform(GENOME.class,width=interval,jump=jump_size,type=2, whole.data=TRUE)
GENOME.class.slide <- F_ST.stats(GENOME.class.slide, mode="nucleotide")
FST_all_slide <- GENOME.class.slide@nuc.F_ST.vs.all
FST_pairwise_slide <- GENOME.class.slide@nuc.F_ST.pairwise

FST_all_slide_d <- as.data.frame(FST_all_slide)
FST_pairwise_slide_d <- as.data.frame(as.vector(FST_pairwise_slide))
#x axis
ids <- length(GENOME.class.slide@region.names)
xaxis <- seq(from = 1, to = ids, by = 1)

#Plot individual populations
for (i in seq_along(population_names))
{
  file_slide <- paste(dir, "_", population_names[i], "_total_FST_per_sliding_window.pdf", sep="")
  slide_plot <- ggplot(FST_all_slide_d, aes(x=xaxis, y=FST_all_slide_d[,i])) + geom_smooth(colour="black", fill="plum") + ggtitle(dir) + xlab("Contig coordinate (kbp)") + ylab("Total FST per interval") + scale_x_continuous(breaks = pretty(xaxis, n = 10))
  ggsave(file_slide, slide_plot)
  #write table with raw data
  slide_table <- paste(dir, "_", population_names[i], "_total_FST_per_sliding_window.txt", sep="")
  write.table(FST_all_slide[,i], file=slide_table, sep="\t",quote=FALSE, col.names=FALSE)
}

#Plot pairwise FST
for (i in seq_along(nrow(FST_pairwise_slide)))
{
  file_hist <- paste(dir, "_", population_names[i], "_pairwise_FST_per_sliding_window", i, ".pdf", sep="")
  slide_plot <- ggplot(FST_pairwise_slide_d, aes(x=xaxis, y=FST_pairwise_slide_d[,i])) + geom_smooth(colour="black", fill="slateblue") + ggtitle(dir) + xlab("Contig coordinate (kbp)") + ylab("Pairwise FST per interval") + scale_x_continuous(breaks = pretty(xaxis, n = 10))
  ggsave(file_slide, slide_plot)
  #write table with raw data
  slide_table <- paste(dir, "_", population_names[i], "_pairwise_FST_per_sliding_window", i, ".txt", sep="")
  fst_table <- cbind(GENOME.class.slide@region.names, FST_pairwise_slide_d[,i])
  write.table(fst_table, file=slide_table, sep="\t",quote=FALSE, col.names=FALSE, row.names=FALSE)
}

}
