#!/bin/bash
input=/home/sobczm/popgen/summary_stats/noA13
vcftools=/home/sobczm/bin/vcftools/bin
#Calculate D, D' and r^2 for SNPs separated by between 1 and 100 kbp
#in non-pathogens
$vcftools/vcftools --vcf Fus2_canu_contigs_unmasked_noA13_filtered.recode.vcf \
--hap-r2 --ld-window-bp-min 1000 --ld-window-bp 100000 \
--indv FOCPG --indv FOCHB6 --indv FOCCB3 --indv FOCA28 --indv FOCD2 --indv FOCA1-2
mv out.hap.ld ld.nonpatho

#start R session
R
library(ggplot2)
table <- "ld.nonpatho"
input <- as.data.frame(read.table(table,header=TRUE))
distance <- input$POS2 - input$POS1
#Plot R^2 against distance
pdf("ld.nonpatho_r2.pdf", width=11, height=8)
qplot(distance, input$R.2)
dev.off()
#Plot D' against distance
pdf("ld.nonpatho_d.pdf", width=11, height=8)
qplot(distance, input$Dprime)
dev.off()
#And histogram
pdf("ld.nonpatho_d_hist.pdf", width=11, height=8)
qplot(input$Dprime, geom="histogram")
dev.off()

#Calculate D, D' and r^2 for SNPs separated by between 1 and 100 kbp
#in pathogens
$vcftools/vcftools --vcf Fus2_canu_contigs_unmasked_noA13_filtered.recode.vcf \
--hap-r2 --ld-window-bp-min 1000 --ld-window-bp 100000 \
--indv FOC55 --indv FOCA23 --indv FOC125 --indv FOCFus2

mv out.hap.ld ld.patho
R
library(ggplot2)
table <- "ld.patho"
input <- as.data.frame(read.table(table,header=TRUE))
distance <- input$POS2 - input$POS1

pdf("ld.patho_r2.pdf", width=11, height=8)
qplot(distance, input$R.2)
dev.off()
#Plot D' against distance
pdf("ld.npatho_d.pdf", width=11, height=8)
qplot(distance, input$Dprime)
dev.off()
#And histogram
pdf("ld.patho_d_hist.pdf", width=11, height=8)
qplot(input$Dprime, geom="histogram")
dev.off()
