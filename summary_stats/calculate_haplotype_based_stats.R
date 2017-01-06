library("PopGenome")
library(ggplot2)

######################## BEFORE RUNNING ################################
#Assign individuals to appropriate populations (or just 1!)
#This script calculates: haplotype-based statistics.
#Output files either per contig or the entire genome (prefix genome_)
##!! Warning, only use on phased genotypes in diploid organisms.
#More than one population needs to be defined, of course!

##!!!!!!!  DIPLOID organisms
##When using diploid organisms and input FASTA files generated using vcf_to_fasta.py, each sample will be artificially
##split into two sequences (<sample_name> + prefix (_1 or _2), for example FOC5_1, FOC5_2), each representing
##one haplotype. Both need to be input below.
nonpatho <- c("FOCA1-2", "FOCA28", "FOCCB3", "FOCD2", "FOCHB6", "FOCPG")
patho <- c("FOCA23", "FOC55", "FOC125", "FOCFus2")
populations <- list(nonpatho, patho)
#Number of populations assigned above.
population_no <- length(populations)
pairs <- choose(population_no,2)
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

  GENOME.class.split <- splitting.data(GENOME.class, subsites="gene")
  GENOME.class.split <- F_ST.stats(GENOME.class.split)
  get.F_ST(GENOME.class.split)

  GENOME.class.slide <- sliding.window.transform(GENOME.class,width=interval,jump=jump_size,type=2, whole.data=TRUE)
  GENOME.class.slide <- F_ST.stats(GENOME.class.slide)
  get.F_ST(GENOME.class.slide)

  ################################################################################
  ###Calculate EXACTLY the same diversity statistiscs, as in calculate_nucleotide_diversity.R and calculate_fst.R
  ###BUT using HAPLOTYPE not NUCLEOTIDE sequences as input.

  Pi <- GENOME.class.split@hap.diversity.within
  Pi_d <- as.data.frame(Pi)

  #Loop over each population: print figure and table with raw data to file
  for (i in seq_along(population_names))
  {
    file_hist <- paste(dir, "_", population_names[i], "_Pi_hap_per_gene.pdf", sep="")
    pi_plot <- ggplot(Pi_d, aes(x=Pi_d[,i])) + geom_histogram(colour="black", fill="blue") + ggtitle(dir) + xlab(expression(paste(pi, " (hap) per gene"))) + ylab("Number of genes") + scale_x_continuous(breaks = pretty(Pi_d[,i], n = 10))
    ggsave(file_hist, pi_plot)
    file_table = paste(dir, "_", population_names[i], "_Pi_hap_per_gene.txt", sep="")
    file_table2 = paste("genome_", population_names[i], "_Pi_hap_per_gene_all.txt", sep="")
    current_gff <- paste(gff, "/", dir, ".gff", sep="")
    gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr=dir, feature=FALSE, extract.gene.names=TRUE)
    Pi_table <- cbind(gene_ids, Pi[,i])
    write.table(Pi_table, file=file_table, sep="\t",quote=FALSE, col.names=FALSE)
    write.table(Pi_table, file=file_table2, sep="\t",quote=FALSE, col.names=FALSE, append=TRUE)

  }

  Pi_persite = GENOME.class.slide@hap.diversity.within
  Pi_persite_d <- as.data.frame(Pi_persite)
  #x axis
  ids <- length(GENOME.class.slide@region.names)
  xaxis <- seq(from = 1, to = ids, by = 1)

  #Plot individual populations
  for (i in seq_along(population_names))
  {
    file_slide <- paste(dir, "_", population_names[i], "_Pi_hap_sliding_window.pdf", sep="")
    slide_plot <- ggplot(Pi_persite_d, aes(x=xaxis, y=Pi_persite_d[,i])) + geom_smooth(colour="black", fill="red") + ggtitle(dir) + xlab("Contig coordinate (kbp)") + ylab(expression(paste(pi, " (hap) per ", interval))) + scale_x_continuous(breaks = pretty(xaxis, n = 10))
    ggsave(file_slide, slide_plot)
    #write table with raw data
    slide_table <- paste(dir, "_", population_names[i], "_Pi_hap_per_sliding_window.txt", sep="")
    write.table(Pi_persite[,i], file=slide_table, sep="\t",quote=FALSE, col.names=FALSE)
  }

  #Plot both populations for comparison
  title <- paste(dir, "Comparison of", population_names[1], "ver.", population_names[2], sep=" ")
  comp_slide_file <- paste(dir, "_Pi_hap_sliding_window_comparison.pdf", sep="")
  slide_comparison <- ggplot(Pi_persite_d, aes(x=xaxis)) + geom_smooth(aes(y=Pi_persite_d[,1]), colour="red") + geom_smooth(aes(y=Pi_persite_d[,2]), colour="blue") + ggtitle(title) + xlab("Contig coordinate (kbp)") + ylab(expression(paste("Average ", pi, " per site"))) + scale_x_continuous(breaks = pretty(xaxis, n = 10))
  ggsave(comp_slide_file, slide_comparison)


  dxy <- GENOME.class.split@hap.diversity.between
  dxy_d <- as.data.frame(as.vector(dxy))
  current_gff <- paste(gff, "/", dir, ".gff", sep="")
  gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr=dir, feature=FALSE, extract.gene.names=TRUE)
  dxy_table <- cbind(GENOME.class.split@region.names, gene_ids, as.vector(dxy))

  #print a histogram of Dxy distribution
  #write table with raw data
  file_hist = paste(dir, "_", "dxy_hap_per_gene.pdf", sep="")
  dxy_plot <- ggplot(dxy_d, aes(x=dxy_d[,1])) + geom_histogram(colour="black", fill="green") + ggtitle(dir) + xlab("Dxy (hap) per gene") + ylab("Number of genes") + scale_x_continuous(breaks = pretty(dxy_d[,1], n = 10))
  ggsave(file_hist, dxy_plot)
  file_table = paste(dir, "_","dxy_hap_per_gene.txt", sep="")
  file_table2 = "genome_dxy_hap_per_gene_all.txt"
  write.table(dxy_table, file=file_table, sep="\t",quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(dxy_table, file=file_table2, sep="\t",quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)

  dxy <- GENOME.class.slide@hap.diversity.between
  dxy_d <- as.data.frame(as.vector(dxy))
  dxy_table <- cbind(GENOME.class.slide@region.names, as.vector(dxy))

  #Plot Dxy across the intervals
  #write table with raw data
  file_slide = paste(dir, "_", "dxy_hap_per_sliding_window.pdf", sep="")
  dxy_plot <- slide_plot <- ggplot(dxy_d, aes(x=xaxis, y=dxy_d[,1])) + geom_smooth(colour="black", fill="green") + ggtitle(dir) + xlab("Contig coordinate (kbp)") + ylab(paste("Dxy (hap) per ", interval, " bp")) + scale_x_continuous(breaks = pretty(xaxis, n = 10))
  ggsave(file_slide, dxy_plot)
  file_table = paste(dir, "_","dxy_hap_per_sliding_window.txt", sep="")
  write.table(dxy_table, file=file_table, sep="\t",quote=FALSE, row.names=FALSE, col.names=FALSE)

  #### Gene-based analysis
  FST_all <- GENOME.class.split@hap.F_ST.vs.all
  FST_all_d <- as.data.frame(FST_all)
  FST_pairwise <- GENOME.class.split@hap.F_ST.pairwise
  FST_pairwise_d <- as.data.frame(as.vector(FST_pairwise))

  for (i in seq_along(population_names))
  {
    file_hist <- paste(dir, "_", population_names[i], "_total_FST_hap_per_gene.pdf", sep="")
    fst_plot <- ggplot(FST_all_d, aes(x=FST_all_d[,i])) + geom_histogram(colour="black", fill="darkseagreen") + ggtitle(dir) + xlab(expression(paste("Total FST (hap) per gene"))) + ylab("Number of genes") + scale_x_continuous(breaks = pretty(FST_all_d[,i], n = 10))
    ggsave(file_hist, fst_plot)
    file_table = paste(dir, "_", population_names[i], "_total_FST_hap_per_gene.txt", sep="")
    file_table2 = paste("genome_", population_names[i], "_total_FST_hap_per_gene_all.txt", sep="")
    current_gff <- paste(gff, "/", dir, ".gff", sep="")
    gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr=dir, feature=FALSE, extract.gene.names=TRUE)
    fst_table <- cbind(gene_ids, FST_all[,i])
    write.table(fst_table, file=file_table, sep="\t",quote=FALSE, col.names=FALSE)
    write.table(fst_table, file=file_table2, sep="\t",quote=FALSE, col.names=FALSE, append=TRUE)
  }

  for (i in seq(pairs))
  {
    file_hist <- paste(dir, "_pairwise_FST_hap_per_gene", i, ".pdf", sep="")
    fst_plot <- ggplot(FST_pairwise_d, aes(x=FST_pairwise_d[,i])) + geom_histogram(colour="black", fill="cadetblue") + ggtitle(dir) + xlab(expression(paste("Pairwise FST per gene"))) + ylab("Number of genes") + scale_x_continuous(breaks = pretty(FST_pairwise_d[,i], n = 10))
    ggsave(file_hist, fst_plot)
    file_table = paste(dir, "_pairwise_FST_hap_per_gene", i, ".txt", sep="")
    file_table2 = paste("genome_pairwise_FST_hap_per_gene_all", i, ".txt", sep="")
    current_gff <- paste(gff, "/", dir, ".gff", sep="")
    gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr=dir, feature=FALSE, extract.gene.names=TRUE)
    fst_table <- cbind(gene_ids, FST_pairwise_d[,i])
    write.table(fst_table, file=file_table, sep="\t",quote=FALSE, col.names=FALSE, row.names=FALSE)
    write.table(fst_table, file=file_table2, sep="\t",quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)
  }

  #### Sliding window analysis (interval)
  FST_all_slide <- GENOME.class.slide@hap.F_ST.vs.all
  FST_pairwise_slide <- GENOME.class.slide@hap.F_ST.pairwise

  FST_all_slide_d <- as.data.frame(FST_all_slide)
  FST_pairwise_slide_d <- as.data.frame(as.vector(FST_pairwise_slide))
  #x axis
  ids <- length(GENOME.class.slide@region.names)
  xaxis <- seq(from = 1, to = ids, by = 1)

  #Plot individual populations
  for (i in seq_along(population_names))
  {
    file_slide <- paste(dir, "_", population_names[i], "_total_FST_hap_per_sliding_window.pdf", sep="")
    slide_plot <- ggplot(FST_all_slide_d, aes(x=xaxis, y=FST_all_slide_d[,i])) + geom_smooth(colour="black", fill="plum") + ggtitle(dir) + xlab("Contig coordinate (kbp)") + ylab("Total FST (hap) per interval") + scale_x_continuous(breaks = pretty(xaxis, n = 10))
    ggsave(file_slide, slide_plot)
    #write table with raw data
    slide_table <- paste(dir, "_", population_names[i], "_total_FST_hap_per_sliding_window.txt", sep="")
    write.table(FST_all_slide[,i], file=slide_table, sep="\t",quote=FALSE, col.names=FALSE)
  }

  #Plot pairwise FST
  for (i in seq(pairs))
  {
    file_hist <- paste(dir, "_pairwise_FST_hap_per_sliding_window", i, ".pdf", sep="")
    slide_plot <- ggplot(FST_pairwise_slide_d, aes(x=xaxis, y=FST_pairwise_slide_d[,i])) + geom_smooth(colour="black", fill="slateblue") + ggtitle(dir) + xlab("Contig coordinate (kbp)") + ylab("Pairwise FST (hap) per interval") + scale_x_continuous(breaks = pretty(xaxis, n = 10))
    ggsave(file_slide, slide_plot)
    #write table with raw data
    slide_table <- paste(dir, "_pairwise_FST_hap_per_sliding_window", i, ".txt", sep="")
    fst_table <- cbind(GENOME.class.slide@region.names, FST_pairwise_slide_d[,i])
    write.table(fst_table, file=slide_table, sep="\t",quote=FALSE, col.names=FALSE, row.names=FALSE)
  }

#FOUR GAMETE test
#per gene
GENOME.class.split <- recomb.stats(GENOME.class.split)
fourgamete_split <- get.recomb(GENOME.class.split)
#per interval
GENOME.class.slide <- recomb.stats(GENOME.class.slide)
fourgamete_slide <- get.recomb(GENOME.class.slide)

#Loop over each population: print figure and table with raw data to file
for (i in seq_along(population_names))
{
  fgt <- unlist(fourgamete_split[i])
  file_hist <- paste(dir, "_", population_names[i], "_4GT_per_gene.pdf", sep="")
  fgt_plot <- ggplot(as.data.frame(fgt), aes(x=as.data.frame(fgt))) + geom_histogram(colour="black", fill="cornsilk") + ggtitle(dir) + xlab("Four gamete test") + ylab("Number of genes") + scale_x_continuous(breaks = pretty(fgt, n = 10))
  ggsave(file_hist, fgt_plot)
  file_table = paste(dir, "_", population_names[i], "_4GT_per_gene.txt", sep="")
  file_table2 = paste("genome_", population_names[i], "_4GT_per_gene_all.txt", sep="")
  current_gff <- paste(gff, "/", dir, ".gff", sep="")
  gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr=dir, feature=FALSE, extract.gene.names=TRUE)
  fgt_table <- cbind(gene_ids, as.data.frame(fourgamete_split[i]))
  write.table(fgt_table, file=file_table, sep="\t",quote=FALSE, col.names=FALSE)
  write.table(fgt_table, file=file_table2, sep="\t",quote=FALSE, col.names=FALSE, append=TRUE)
}

for (i in seq_along(population_names))
{
  fgt <- unlist(fourgamete_slide[i])
  file_slide <- paste(dir, "_", population_names[i], "_4GT_per_sliding_window.pdf", sep="")
  slide_plot <- ggplot(as.data.frame(fgt), aes(x=xaxis, y=as.data.frame(fgt))) + geom_smooth(colour="black", fill="greenyellow") + ggtitle(dir) + xlab("Contig coordinate (kbp)") + ylab("Four gamete test") + scale_x_continuous(breaks = pretty(xaxis, n = 10))
  ggsave(file_slide, slide_plot)
  #write table with raw data
  slide_table <- paste(dir, "_", population_names[i], "_4GT_per_sliding_window.txt", sep="")
  write.table(as.data.frame(fourgamete_slide[i]), file=slide_table, sep="\t",quote=FALSE, col.names=FALSE)
}

}

for (i in seq_along(population_names))
{
  #Pi (hap)
  file_table2 = paste("genome_", population_names[i], "_Pi_hap_per_gene_all.txt", sep="")
  x <- as.data.frame(read.delim(file_table2))
  file_hist <- paste("genome_", population_names[i], "_Pi_hap_per_gene_all.pdf", sep="")
  pi_plot <- ggplot(x, aes(x=x[,3])) + geom_histogram(colour="black", fill="blue") + xlab(expression(paste(pi, " (hap) per gene"))) + ylab("Number of genes") + scale_x_continuous(breaks = pretty(x[,3], n = 10))
  ggsave(file_hist, pi_plot)
  #Total FST (hap)
  file_table2 = paste("genome_", population_names[i], "_total_FST_hap_per_gene_all.txt", sep="")
  x <- as.data.frame(read.delim(file_table2))
  file_hist <- paste("genome_", population_names[i], "_total_FST_hap_per_gene_all.pdf", sep="")
  fst_plot <- ggplot(x, aes(x=x[,3])) + geom_histogram(colour="black", fill="darkseagreen") + xlab(expression(paste("Total FST (hap) per gene"))) + ylab("Number of genes") + scale_x_continuous(breaks = pretty(x[,3], n = 10))
  ggsave(file_hist, fst_plot)
  #Four gamete test
  file_table2 = paste("genome_", population_names[i], "_4GT_per_gene_all.txt", sep="")
  x <- as.data.frame(read.delim(file_table2))
  file_hist <- paste("genome_", population_names[i], "_4GT_per_gene_all.pdf", sep="")
  fgt_plot <- ggplot(x, aes(x=x[,3])) + geom_histogram(colour="black", fill="cornsilk") + xlab("Four gamete test") + ylab("Number of genes") + scale_x_continuous(breaks = pretty(x[,3], n = 10))
  ggsave(file_hist, fgt_plot)
}

for (i in seq(pairs))
{
#Pairwise FST (hap)
file_table2 = paste("genome_pairwise_FST_hap_per_gene_all", i, ".txt", sep="")
x <- as.data.frame(read.delim(file_table2))
file_hist <- paste("genome_pairwise_FST_hap_per_gene_all", i, ".pdf", sep="")
fst_plot <- ggplot(x, aes(x=x[,2])) + geom_histogram(colour="black", fill="cadetblue") + xlab(expression(paste("Pairwise FST (hap) per gene"))) + ylab("Number of genes") + scale_x_continuous(breaks = pretty(x[,2], n = 10))
ggsave(file_hist, fst_plot)
}

#Dxy (hap)
file_table2 = "genome_dxy_hap_per_gene_all.txt"
x <- as.data.frame(read.delim(file_table2))
file_hist = "genome_dxy_hap_per_gene_all.pdf"
dxy_plot <- ggplot(x, aes(x=x[,3])) + geom_histogram(colour="black", fill="green") + xlab("Dxy (hap) per gene") + ylab("Number of genes") + scale_x_continuous(breaks = pretty(x[,3], n = 10))
ggsave(file_hist, dxy_plot)

#Linkage disequilibrium stats
#GENOME.class.slide <- calc.R2(GENOME.class.slide)
#GENOME.class.split<- calc.R2(GENOME.class.split)
#GENOME.class.slide@region.stats@linkage.disequilibrium
#GENOME.class.split@region.stats@linkage.disequilibrium

#Loop over each population: print table with raw data to file
#for (i in seq_along(population_names))
#{
#for (j in seq_along(gene_ids))
#{
#  current_locus <-
#  current_locus_l <- length(GENOME.class.split@region.stats@linkage.disequilibrium[[j]][[i]]) / 3
#  for (z in range(1, current_locus_l))
#  {
#  if (!is.null(current_locus))
#  {
#  if (!is.na(current_locus))
#  recomb_table <- cbind(gene_ids[j], current_locus[1,z], current_locus[2,z], current_locus[3,z])
#  file_table = paste(dir, "_", population_names[i], "_r2.txt", sep="")
#  write.table(recomb_table, file=file_table, sep="\t",quote=FALSE, col.names=FALSE, append=TRUE)
#  }
#  }
#  }
#}

#[[x]][[y]] x:region, y:population
