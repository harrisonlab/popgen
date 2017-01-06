library("PopGenome")
library(ggplot2)
######################## BEFORE RUNNING ################################
#Assign individuals to appropriate populations (or just 1!)
#However, if using just 1 populations, interpopulation statistic (Dxy)
#of genetic variation cannot be calculated, and need to carry out only analyses
#A, B  but not C, D.
#When all coding sites input, use E to calculate Pi(nonsyn)/Pi(syn) sites.
#Output files either per contig or the entire genome (prefix genome_)
##!!!!!!!  DIPLOID organisms
##When using diploid organisms and input FASTA files generated using vcf_to_fasta.py, each sample will be artificially
##split into two sequences (<sample_name> + prefix (_1 or _2), for example FOC5_1, FOC5_2), each representing
##one haplotype. Both need to be input below.
nonpatho <- c("FOCA1-2", "FOCA28", "FOCCB3", "FOCD2", "FOCHB6", "FOCPG")
patho <- c("FOCA23", "FOC55", "FOC125", "FOCFus2")
populations <- list(nonpatho, patho)
#Number of populations assigned above.
population_no <- length(populations)
population_names <- c("nonpatho", "patho")
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

#########################################################
#A) calculate Pi (Nei, 1987) for all sites in a given gene.

#Split by and retain only genes for the analysis
GENOME.class.split <- splitting.data(GENOME.class, subsites="gene")
GENOME.class.split <- diversity.stats(GENOME.class.split, pi=TRUE)

#Divide Pi per number of sites in the gene to calculate value per site
Pi <- GENOME.class.split@Pi / GENOME.class.split@n.sites
Pi_d <- as.data.frame(Pi)

#Loop over each population: print figure and table with raw data to file
for (i in seq_along(population_names))
{
file_hist <- paste(dir, "_", population_names[i], "_Pi_per_gene.pdf", sep="")
pi_plot <- ggplot(Pi_d, aes(x=Pi_d[,i])) + geom_histogram(colour="black", fill="blue") + ggtitle(dir) + xlab(expression(paste("Average ", pi, " per site"))) + ylab("Number of genes") + scale_x_continuous(breaks = pretty(Pi_d[,i], n = 10))
ggsave(file_hist, pi_plot)
file_table = paste(dir, "_", population_names[i], "_Pi_per_gene.txt", sep="")
file_table2 = paste("genome_", population_names[i], "_Pi_per_gene_all.txt", sep="")
current_gff <- paste(gff, "/", dir, ".gff", sep="")
gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr=dir, feature=FALSE, extract.gene.names=TRUE)
Pi_table <- cbind(gene_ids, Pi[,i])
write.table(Pi_table, file=file_table, sep="\t",quote=FALSE, col.names=FALSE)
#Table with genome-wide results:
write.table(Pi_table, file=file_table2, sep="\t",quote=FALSE, col.names=FALSE, append=TRUE)
}
###########################################################
#B) calculate Pi (Nei, 1987) in a sliding-window for all sites over a given interval
GENOME.class.slide <- sliding.window.transform(GENOME.class,width=interval,jump=jump_size,type=2, whole.data=TRUE)
GENOME.class.slide <- diversity.stats(GENOME.class.slide, pi=TRUE)
#plot the results for all populations over one figure
Pi_persite = GENOME.class.slide@Pi / interval
Pi_persite_d <- as.data.frame(Pi_persite)
#x axis
ids <- length(GENOME.class.slide@region.names)
xaxis <- seq(from = 1, to = ids, by = 1)

#Plot individual populations
for (i in seq_along(population_names))
{
file_slide <- paste(dir, "_", population_names[i], "_Pi_sliding_window.pdf", sep="")
slide_plot <- ggplot(Pi_persite_d, aes(x=xaxis, y=Pi_persite_d[,i])) + geom_smooth(colour="black", fill="red") + ggtitle(dir) + xlab("Contig coordinate (kbp)") + ylab(expression(paste("Average ", pi, " per site"))) + scale_x_continuous(breaks = pretty(xaxis, n = 10))
ggsave(file_slide, slide_plot)
#write table with raw data
slide_table <- paste(dir, "_", population_names[i], "_Pi_per_sliding_window.txt", sep="")
write.table(Pi_persite[,i], file=slide_table, sep="\t",quote=FALSE, col.names=FALSE)
}

#Plot both populations for comparison
title <- paste(dir, "Comparison of", population_names[1], "ver.", population_names[2], sep=" ")
comp_slide_file <- paste(dir, "_Pi_sliding_window_comparison.pdf", sep="")
slide_comparison <- ggplot(Pi_persite_d, aes(x=xaxis)) + geom_smooth(aes(y=Pi_persite_d[,1]), colour="red") + geom_smooth(aes(y=Pi_persite_d[,2]), colour="blue") + ggtitle(title) + xlab("Contig coordinate (kbp)") + ylab(expression(paste("Average ", pi, " per site"))) + scale_x_continuous(breaks = pretty(xaxis, n = 10))
ggsave(comp_slide_file, slide_comparison)

############################################################
#C) Calculate Dxy for all sites in a given gene, if more than 1 population analysed.
#All possible pairwise contrasts
GENOME.class.split <- F_ST.stats(GENOME.class.split)
FST_results <- get.F_ST(GENOME.class.split)
dxy <- GENOME.class.split@nuc.diversity.between / GENOME.class.split@n.sites
dxy_d <- as.data.frame(as.vector(dxy))
current_gff <- paste(gff, "/", dir, ".gff", sep="")
gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr=dir, feature=FALSE, extract.gene.names=TRUE)
dxy_table <- cbind(GENOME.class.split@region.names, gene_ids, as.vector(dxy))

#print a histogram of Dxy distribution
#write table with raw data
file_hist = paste(dir, "_", "dxy_per_gene.pdf", sep="")
dxy_plot <- ggplot(dxy_d, aes(x=dxy_d[,1])) + geom_histogram(colour="black", fill="green") + ggtitle(dir) + xlab("Average Dxy per gene") + ylab("Number of genes") + scale_x_continuous(breaks = pretty(dxy_d[,1], n = 10))
ggsave(file_hist, dxy_plot)
file_table = paste(dir, "_","dxy_per_gene.txt", sep="")
file_table2 = "genome_dxy_per_gene_all.txt"
write.table(dxy_table, file=file_table, sep="\t",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(dxy_table, file=file_table2, sep="\t",quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
############################################################
#D) Calculate Dxy for all sites in a sliding window analysis, if more than 1 population analysed.
#All possible pairwise contrasts
GENOME.class.slide <- F_ST.stats(GENOME.class.slide)
FST_results <- get.F_ST(GENOME.class.slide)
dxy <- GENOME.class.slide@nuc.diversity.between / GENOME.class.slide@n.sites
dxy_d <- as.data.frame(as.vector(dxy))
dxy_table <- cbind(GENOME.class.slide@region.names, as.vector(dxy))

#Plot Dxy across the intervals
#write table with raw data
file_slide = paste(dir, "_", "dxy_per_sliding_window.pdf", sep="")
dxy_plot <- slide_plot <- ggplot(dxy_d, aes(x=xaxis, y=dxy_d[,1])) + geom_smooth(colour="black", fill="green") + ggtitle(dir) + xlab("Contig coordinate (kbp)") + ylab(paste("Average Dxy per ", interval, " bp")) + scale_x_continuous(breaks = pretty(xaxis, n = 10))
ggsave(file_slide, dxy_plot)
file_table = paste(dir, "_","dxy_per_sliding_window.txt", sep="")
write.table(dxy_table, file=file_table, sep="\t",quote=FALSE, row.names=FALSE, col.names=FALSE)

##############################################################
#E) When a dataset containing all types of coding sites loaded, calculate Pi(nonsyn)/Pi(syn) over each
# gene and over a given interval in the genome.
#Gene-based
GENOME.class.split.nonsyn <- diversity.stats(GENOME.class.split, pi=TRUE, subsites="nonsyn")
GENOME.class.split.syn <- diversity.stats(GENOME.class.split, pi=TRUE, subsites="syn")
#Interval-based
GENOME.class.slide.nonsyn <- diversity.stats(GENOME.class.slide, pi=TRUE, subsites="nonsyn")
GENOME.class.slide.syn <- diversity.stats(GENOME.class.slide, pi=TRUE, subsites="syn")

## Print output (gene-based)
#Plot individual populations (gene-based)
#Divide Pi per number of sites in the gene to calculate value per site
Pi_ns <- GENOME.class.split.nonsyn@Pi / GENOME.class.split.syn@Pi / GENOME.class.split@n.sites
Pi_ns_d <- as.data.frame(Pi_ns)

#Loop over each population: print figure and table with raw data to file
for (i in seq_along(population_names))
{
  #Check in case all values 0
  Pi_ns_len <- length(Pi_ns[,i])
  Pi_ns_len_na <- sum(sapply(Pi_ns[,i], is.na))
  if (Pi_ns_len > Pi_ns_len_na){
  file_hist <- paste(dir, "_", population_names[i], "_Pi_n_s_per_gene.pdf", sep="")
  pi_plot <- ggplot(Pi_ns_d, aes(x=Pi_ns_d[,i])) + geom_histogram(colour="black", fill="coral") + ggtitle(dir) + xlab(expression(paste("Average ", pi, "ns/", pi, "s", " per site"))) + ylab("Number of genes") + scale_x_continuous(breaks = pretty(Pi_ns_d[,i], n = 10))
  ggsave(file_hist, pi_plot)}
  file_table = paste(dir, "_", population_names[i], "_Pi_n_s_per_gene.txt", sep="")
  file_table2 = paste("genome_", population_names[i], "_Pi_n_s_per_gene_all.txt", sep="")
  current_gff <- paste(gff, "/", dir, ".gff", sep="")
  gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr=dir, feature=FALSE, extract.gene.names=TRUE)
  Pi_table <- cbind(gene_ids, Pi_ns[,i])
  write.table(Pi_table, file=file_table, sep="\t",quote=FALSE, col.names=FALSE)
  write.table(Pi_table, file=file_table2, sep="\t",quote=FALSE, col.names=FALSE, append=TRUE)
}

## Print output (interval-based)
#plot the results for all populations over one figure
Pi_ns_persite <- GENOME.class.slide.nonsyn@Pi / GENOME.class.slide.syn@Pi / interval
Pi_ns_persite_d <- as.data.frame(Pi_ns_persite)
#x axis
ids <- length(GENOME.class.slide@region.names)
xaxis <- seq(from = 1, to = ids, by = 1)

#Plot individual populations
for (i in seq_along(population_names))
{
  file_slide <- paste(dir, "_", population_names[i], "_Pi_n_s_sliding_window.pdf", sep="")
  Pi_ns_len <- length(Pi_ns_persite[,i])
  Pi_ns_len_na <- sum(sapply(Pi_ns_persite[,i], is.na))
  if (Pi_ns_len > Pi_ns_len_na){
  slide_plot <- ggplot(Pi_ns_persite_d, aes(x=xaxis, y=Pi_ns_persite_d[,i])) + geom_smooth(colour="black", fill="darkviolet") + ggtitle(dir) + xlab("Contig coordinate (kbp)") + ylab(expression(paste("Average ", pi, "ns/", pi, "s", " per site"))) + scale_x_continuous(breaks = pretty(xaxis, n = 10))
  ggsave(file_slide, slide_plot)}
  #write table with raw data
  slide_table <- paste(dir, "_", population_names[i], "_Pi_n_s_per_sliding_window.txt", sep="")
  write.table(Pi_ns_persite[,i], file=slide_table, sep="\t",quote=FALSE, col.names=FALSE)
}

#Plot both populations for comparison
title <- paste(dir, "Comparison of", population_names[1], "ver.", population_names[2], sep=" ")
comp_slide_file <- paste(dir, "_Pi_n_s_sliding_window_comparison.pdf", sep="")
if (Pi_ns_len > Pi_ns_len_na){
slide_comparison <- ggplot(Pi_ns_persite_d, aes(x=xaxis)) + geom_smooth(aes(y=Pi_ns_persite_d[,1]), colour="deeppink") + geom_smooth(aes(y=Pi_ns_persite_d[,2]), colour="lightskyblue") + ggtitle(title) + xlab("Contig coordinate (kbp)") + ylab(expression(paste("Average ", pi, "ns/", pi, "s", " per site"))) + scale_x_continuous(breaks = pretty(xaxis, n = 10))
ggsave(comp_slide_file, slide_comparison)}

}

###Plot genome-wide histograms
for (i in seq_along(population_names))
{
#Pi table
file_table2 <- paste("genome_", population_names[i], "_Pi_per_gene_all.txt", sep="")
x <- as.data.frame(read.delim(file_table2))
file_hist <- paste("genome_", population_names[i], "_Pi_per_gene_all.pdf", sep="")
pi_plot <- ggplot(x, aes(x=x[,3])) + geom_histogram(colour="black", fill="blue") + xlab(expression(paste("Average ", pi, " per site"))) + ylab("Number of genes") + scale_x_continuous(breaks = pretty(x[,3], n = 10))
ggsave(file_hist, pi_plot)

#Pi nonsyn/syn
file_table2 = paste("genome_", population_names[i], "_Pi_n_s_per_gene_all.txt", sep="")
x <- as.data.frame(read.delim(file_table2))
file_hist <- paste("genome_", population_names[i], "_Pi_n_s_per_gene_all.pdf", sep="")
pi_plot <- ggplot(x, aes(x=x[,3])) + geom_histogram(colour="black", fill="coral") + xlab(expression(paste("Average ", pi, "ns/", pi, "s", " per site"))) + ylab("Number of genes") + scale_x_continuous(breaks = pretty(x[,3], n = 10))
ggsave(file_hist, pi_plot)
}

#Dxy table
file_table2 = "genome_dxy_per_gene_all.txt"
x <- as.data.frame(read.delim(file_table2))
file_hist = "genome_dxy_per_gene_all.pdf"
dxy_plot <- ggplot(x, aes(x=x[,3])) + geom_histogram(colour="black", fill="green") + xlab("Average Dxy per gene") + ylab("Number of genes") + scale_x_continuous(breaks = pretty(x[,3], n = 10))
ggsave(file_hist, dxy_plot)
