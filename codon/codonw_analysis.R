###Inputs
# x - codonw output file with all genes
# files - a list of files with specific gene subsets from codonw output contained in x

x <- "125_final_genes_combined.cdna_pass_one_notrans_only_an.out"

files <- c("125_CAZY_gene_table_125_final_genes_combined.cdna_pass_one_notrans_only_an",
           "125_CAZY_secreted_gene_table_125_final_genes_combined.cdna_pass_one_notrans_only_an",
           "125_genes_in_2kb_mimp_gene_table_125_final_genes_combined.cdna_pass_one_notrans_only_an",
           "125_genes_in_2kb_mimp_secreted_gene_table_125_final_genes_combined.cdna_pass_one_notrans_only_an",
           "F.oxysporum_fsp_cepae_125_EffectorP_secreted_gene_table_125_final_genes_combined.cdna_pass_one_notrans_only_an"
            )

output = read.table(x, header=T, row.names = 1, stringsAsFactors=F)
#CAI = (output[, 5])
#hist(CAI)
#FOP = (output[, 7])
#hist(FOP)

pdf("CAI_histogram.PDF", width=11, height=8)
hist(output$CAI, breaks=50, col="lightgreen")
dev.off()
pdf("FOP_histogram.PDF", width=11, height=8)
hist(output$Fop, breaks=50, col="pink")
dev.off()

Fop_n <- as.numeric(as.character(output$Fop))
GC_n <- as.numeric(as.character(output$GC))
GC3s_n <- as.numeric(as.character(output$GC3s))
Nc_n <- as.numeric(as.character(output$Nc))
CAI_n <- as.numeric(as.character(output$CAI))
CBI_n <- as.numeric(as.character(output$CBI))

#Read in the Correspondence Analysis results
corr = read.table("genes.coa", header=T, row.names = 1, stringsAsFactors=F)
#Axis 1 vs Nc
prin <- as.numeric(as.character(corr$Axis1))
Nc <- as.numeric(as.character(output$Nc))

pdf("scatterplots.PDF", width=11, height=8)
plot(output$GC3s, output$Nc, xlab = "GC3s",
     ylab = "Nc", main = "Correspondence analysis", cex = 0.1)

plot(output$GC3s, output$Fop, xlab = "Fop",
     ylab = "GC3s", main = "Correspondence analysis", cex = 0.1)
abline(lm(output$Fop ~ output$GC3s))

plot(output$GC, output$Nc, xlab = "Nc",
     ylab = "GC", main = "Correspondence analysis", cex = 0.1)

plot(output$GC, output$Fop, xlab = "Fop",
     ylab = "GC", main = "Correspondence analysis", cex = 0.1)
abline(lm(output$Fop ~ output$GC))

plot(corr$Axis1, corr$Axis2, xlab = "Axis1",
     ylab = "Axis2", main = "Correspondence analysis", cex = 0.1)

plot(corr$Axis1, output$Nc, xlab = "Axis1",
     ylab = "Nc", main = "Correspondence analysis", cex = 0.1)

plot(corr$Axis1, output$Fop, xlab = "Axis1",
     ylab = "Fop", main = "Correspondence analysis", cex = 0.1)
abline(lm(output$Fop ~ corr$Axis1))
dev.off()

sink(file = "correlation.txt", append = T, type = c("output", "message"), split = F)
cor.test(Fop_n, GC_n, method="spearman")
cor.test(Fop_n, GC3s_n, method="spearman")
cor.test(Nc_n, GC_n, method="spearman")
cor.test(Nc_n, GC3s_n, method="spearman")
cor.test(corr$Axis1, Nc, method="spearman")
sink()


#Print summary stats for all genes
data_frame <- c("Name", "n", "CAI", "", "CBI", "", "FOP", "", "Nc", "", "GC3s","", "GC", "")
srow <- c("", "", "mean", "p-value", "mean", "p-value", "mean", "p-value", "mean", "p-value", "mean", "p-value", "mean", "p-value")
data_frame <- rbind(data_frame, srow)
all <- c(x, length(CAI_n), mean(CAI_n), NA, mean(CBI_n), NA, mean(Fop_n), NA, mean(Nc_n, na.rm=TRUE), NA, mean(GC3s_n), NA, mean(GC_n), NA)
data_frame <- rbind(data_frame, all)
###Compare values of individual gene sets with all genes with t-test and print summary stats
for (file in files)
{
  effectors <- read.table(file, header=F, row.names = 1, stringsAsFactors=F)
#Compare CAI
cai_t <- t.test(CAI_n, as.numeric(as.character(effectors[,5])))
#Compare CBI
cbi_t <- t.test(CBI_n, as.numeric(as.character(effectors[,6])))
#compare Fop
fop_t <- t.test(Fop_n, as.numeric(as.character(effectors[,7])))
#Compare Nc
Nc_e <- as.numeric(as.character(effectors[,8]))
Nc_t <- t.test(Nc_n, Nc_e)

#Compare GC3s
gc3s_t <- t.test(GC3s_n, as.numeric(as.character(effectors[,9])))
#Compare GC
gc_t <- t.test(GC_n, as.numeric(as.character(effectors[,10])))
eff <- c(file, length(effectors[,5]), mean(effectors[,5]), cai_t$p.value, mean(effectors[,6]), cbi_t$p.value, mean(effectors[,7]), fop_t$p.value, mean(Nc_e, na.rm=TRUE), Nc_t$p.value, mean(effectors[,9]), gc3s_t$p.value, mean(effectors[,10]), gc_t$p.value)
data_frame <- rbind(data_frame, eff)
}
write.table(data_frame, file="codonw_summary.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
