#!/bin/bash
scripts=/home/sobczm/bin/popgen/summary_stats
input=/home/sobczm/popgen/summary_stats
#################### 1) Filter SNPs to retain only biallelic SNPs. Furthermore, keep only SNPs with max 5% missing genotypes (but can change that with --max-missing option below.
#!! It may be necessary to also pre-filter for samples which were poorly sequenced/
#aligned beforehand to avoid removal of too many potentially informative SNPs.
for filename in $input/outliers/*/*.vcf
do
vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf $filename --max-missing 0.95 --mac 1 --min-alleles 2 --max-alleles 2 --recode --out ${filename%.vcf}_bi_filtered
done

###Warning! Following the inspection of scree plot after the first run, 
#the R script then requires manual changing of the final_k parameter 
#which is the k value best explaining the population structure.

input_hap=$input/outliers/vinequalis/Ash_farm_172_pacbio_contigs_unmasked_no118_bi_filtered.recode.vcf
input_dip=$input/outliers/pfrag/95m_contigs_unmasked_UK123_bi_filtered.recode.vcf


qsub $scripts/sub_pcadapt.sh $input_hap 1
qsub $scripts/sub_pcadapt.sh $input_dip 2
