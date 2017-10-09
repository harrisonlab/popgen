#!/bin/bash
scripts=/home/sobczm/bin/popgen/snp
input=/home/sobczm/popgen/snp/snp_chip/everbearer_gwas
cd $input

#Extract all possible samples, including duplicates.
python $scripts/ananassa_genotypes_db.py all_everbearer_samples.txt all_everbearer_samples.out
python $scripts/ananassa_genotypes_vcf.py all_everbearer_samples.out istraw90_vesca_v2.0_snp_positions.gff3
#How many SNPs left if only 5% missing data allowed?
$vcftools/vcftools --vcf all_everbearer_samples.out.vcf --max-missing 0.95 --recode --out all_everbearer_sample_max005
#kept 16419 out of a possible 95331 Sites
#Sort the resulting VCF files and fix so that compatible with Plink.
vcftools=/home/sobczm/bin/vcftools/bin
for a in *.vcf
do
$vcftools/vcf-sort $a | sed 's/LG//' | sed 's/Unknown/0/' >${a%.vcf}_sorted.vcf
done
#Calculate kinship in PLINK.
plink --vcf all_everbearer_sample_max005.recode_sorted.vcf --genome --out all_everbearer_kinship
cat all_everbearer_kinship.genome | awk '{$1=$1;print}' OFS='\t' >temp
mv temp all_everbearer_kinship.genome
#Plot my distance matrix
$scripts/similarity_percentage.py all_everbearer_sample_max005.recode_sorted.vcf
#Visualise the output as heatmap and clustering dendrogram
Rscript --vanilla $scripts/distance_matrix.R all_everbearer_sample_max005.recode_sorted_distance.log