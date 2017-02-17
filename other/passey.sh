#!/bin/bash
scripts=/home/sobczm/bin/popgen/snp
input=/home/sobczm/popgen/other/passey

#Calculate NJ tree based on SNPs.

#First eliminate sample 118 as it's mostly missing data - poorly sequenced

vcflib=/home/sobczm/bin/vcflib/bin
$vcflib/vcfremovesamples Ash_farm_172_pacbio_contigs_unmasked.vcf 118 >Ash_farm_172_pacbio_contigs_unmasked_no118.vcf

#Now filter the SNPs to retain those with no missing data
#in any individual.
qsub $scripts/sub_vcf_parser.sh Ash_farm_172_pacbio_contigs_unmasked_no118.vcf 40 30 10 30 1 Y

#Prepare an NJ tree
$scripts/nj_tree.sh Ash_farm_172_pacbio_contigs_unmasked_no118_filtered.vcf 1