#!/bin/bash
input=/home/sobczm/popgen/snp/snp_calling/multithreaded
scripts=/home/sobczm/bin/popgen/snp

#Only retain biallelic high-quality SNPS with no missing data for genetic analyses.
cd $input
$scripts/vcf_parser_haploid.py --i Fus2_canu_contigs_unmasked.vcf

#Calculate the index for percentage of shared SNP alleles between the individs.
$scripts/similarity_percentage.py Fus2_canu_contigs_unmasked_filtered.vcf
