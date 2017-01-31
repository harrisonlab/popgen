library(PopGenome)
library(ggplot2)

######################## BEFORE RUNNING ################################
#Assign individuals to appropriate populations (or just 1!)
#This script calculates: Fay Wu's H and requires an outgroup species (coded as 'ancestral' population below)
#Output files either per contig or the entire genome (prefix genome_)
##!!!!!!!  DIPLOID organisms
##When using diploid organisms and input FASTA files generated using vcf_to_fasta.py, each sample will be artificially
##split into two sequences (<sample_name> + prefix (_1 or _2), for example FOC5_1, FOC5_2), each representing
##one haplotype. Both need to be input below.
Pfrag <- c("Bc1_1", "Bc1_2", "Nov5_1", "Nov5_2", "Bc16_1", "Bc16_2", "A4_1", "A4_2", "Nov27_1", "Nov27_2", "Nov9_1", "Nov9_2", "Nov71_1", "Nov71_2")
Prubi <- c("SCRP324_1", "SCRP324_2", "SCRP333_1", "SCRP333_2")
#Assign outgroup samples to the "ancestral" population. The population name "ancestral" should
#not be changed as it evoked below on line 35.
ancestral <- c("ancestral_1_1", "ancestral_1_2", "ancestral_2_1", "ancestral_2_2")
populations <- list(Pfrag, Prubi)
#Number of populations assigned above.
population_no <- length(populations)
population_names <- c("Pfrag", "Prubi")
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
  GENOME.class <- set.outgroup(GENOME.class, new.outgroup=ancestral)

  #Calculate neutrality stats over genes
  GENOME.class.split <- splitting.data(GENOME.class, subsites="gene")
  GENOME.class.split <- neutrality.stats(GENOME.class.split)
  FayWuH <- GENOME.class.split@Fay.Wu.H
  FayWuH_d <- as.data.frame(FayWuH)
  for (i in seq_along(population_names))
  {
    file_hist <- paste(dir, "_", population_names[i], "_FayWuH_per_gene.pdf", sep="")
    FayWuH_len <- length(FayWuH_d[,i])
    FayWuH_na <- sum(sapply(FayWuH[,i], is.na))
    if (FayWuH_len > FayWuH_na){
      faywuh_plot <- ggplot(FayWuH_d, aes(x=FayWuH_d[,i])) + geom_histogram(colour="black", fill="thistle4") + ggtitle(dir) + xlab("Fay & Wu's H") + ylab("Number of genes") + scale_x_continuous(breaks = pretty(FayWuH_d[,i], n = 10))
      ggsave(file_hist, faywuh_plot)}
    file_table = paste(dir, "_", population_names[i], "_FayWuH_per_gene.txt", sep="")
    file_table2 = paste("genome_", population_names[i], "_FayWuH_per_gene_all.txt", sep="")
    current_gff <- paste(gff, "/", dir, ".gff", sep="")
    gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr=dir, feature=FALSE, extract.gene.names=TRUE)
    faywuh_table <- cbind(gene_ids, FayWuH[,i])
    write.table(faywuh_table, file=file_table, sep="\t",quote=FALSE, col.names=FALSE)
    write.table(faywuh_table, file=file_table2, sep="\t",quote=FALSE, col.names=FALSE, append=TRUE)
  }
  
}

###Plot genome-wide histograms
for (i in seq_along(population_names))
{
  #Tajima's D table
  file_table2 <- paste("genome_", population_names[i], "_FayWuH_per_gene_all.txt", sep="")
  x <- as.data.frame(read.delim(file_table2))
  file_hist <- paste("genome_", population_names[i], "_FayWuH_per_gene_all.pdf", sep="")
  faywuh_plot <- ggplot(x, aes(x=x[,3])) + geom_histogram(colour="black", fill="thistle1") + xlab("Fay & Wu's H") + ylab("Number of genes") + scale_x_continuous(breaks = pretty(x[,3], n = 10))
  ggsave(file_hist, faywuh_plot)
}
