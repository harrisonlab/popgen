#!/bin/bash
input=/home/sobczm/popgen/summary_stats/noA13
scripts=/home/sobczm/bin/popgen/summary_stats
vcftools=/home/sobczm/bin/vcftools/bin

#Calculate D, D' and r^2 for SNPs separated by between 1 and 100 kbp
#in non-pathogens (program calculates the stats using only the individuals
#listed after "--indv" switch)
#In order to calculate all versus all SNP comparison, remove options --ld-window-bp-min
#and --ld-window-bp
$vcftools/vcftools --vcf Fus2_canu_contigs_unmasked_noA13_filtered.recode.vcf \
--hap-r2 --ld-window-bp-min 1000 --ld-window-bp 100000 \
--indv FOCPG --indv FOCHB6 --indv FOCCB3 --indv FOCA28 --indv FOCD2 --indv FOCA1-2
mv out.hap.ld ld.nonpatho

qsub $scripts/sub_plot_ld.sh ld.nonpatho

#Calculate D, D' and r^2 for SNPs separated by between 1 and 100 kbp
#in pathogens (program calculates the stats using only the individuals
#listed after "--indv" switch)

$vcftools/vcftools --vcf Fus2_canu_contigs_unmasked_noA13_filtered.recode.vcf \
--hap-r2 --ld-window-bp-min 1000 --ld-window-bp 100000 \
--indv FOC55 --indv FOCA23 --indv FOC125 --indv FOCFus2

mv out.hap.ld ld.patho

#Plot D' and r2 versus SNP physical distance, histogram of D' values
qsub $scripts/sub_plot_ld.sh ld.patho
