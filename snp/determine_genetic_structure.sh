#!/bin/bash
input=/home/sobczm/popgen/snp/snp_calling/multithreaded
scripts=/home/sobczm/bin/popgen/snp

#Only retain biallelic high-quality SNPS with no missing data for genetic analyses.
cd $input
qsub $scripts/sub_vcf_parser.sh Fus2_canu_contigs_unmasked.vcf

#In some organisms, may want to thin (subsample) SNPs in high linkage diseqilibrium down to
#1 SNP  per e.g. 10 kbp just for the population structure analyses.
#vcftools=/home/sobczm/bin/vcftools/bin
#$vcftools/vcftools --vcf $input_vcf --thin 10000 --recode --out ${input_vcf%.vcf}_thinned

#General VCF stats (remember that vcftools needs to have the PERL library exported)
perl /home/sobczm/bin/vcftools/bin/vcf-stats \
Fus2_canu_contigs_unmasked.vcf >Fus2_canu_contigs_unmasked.stat
perl /home/sobczm/bin/vcftools/bin/vcf-stats \
Fus2_canu_contigs_unmasked_filtered.vcf >Fus2_canu_contigs_unmasked_filtered.stat
#Calculate the index for percentage of shared SNP alleles between the individs.
$scripts/similarity_percentage.py Fus2_canu_contigs_unmasked_filtered.vcf
#Visualise the output as heatmap and clustering dendrogram
Rscript --vanilla $scripts/distance_matrix.R Fus2_canu_contigs_unmasked_filtered_distance.log
#Carry out PCA and plot the results
Rscript --vanilla $scripts/pca.R Fus2_canu_contigs_unmasked_filtered.vcf
#Calculate an NJ tree based on all the SNPs. Outputs a basic diplay of the tree, plus a Newick file to be used
#for displaying the tree in FigTree and beautifying it.
#First argument - input VCF file, second argument - ploidy of the samples
$scripts/nj_tree.sh Fus2_canu_contigs_unmasked_filtered.vcf 1
#DAPC and AMOVA analysis
Rscript --vanilla $popgen/snp/amova_dapc.R
