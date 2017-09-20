#!/bin/bash
scripts=/home/sobczm/bin/popgen/snp
input=/home/sobczm/popgen/snp/snp_chip/crown_rot_gwas
cd $input

#First, running plink and TASSEL ithout imputation of genotypes - as we don't really know if we can impute them correctly, anyway - same problem as phasing.

#Copy over the VCF file with target individuals (retaining a couple of EMxFE individuals, but no sample duplicates).
cp /home/sobczm/popgen/snp/snp_chip/crown_rot_local/QTL_positions/sample_ids_crown_rot.out_nodup.vcf $input
#Convert VCF to plink and retain only informative SNPs with MAF of at least 0.05. 
#First, plink does not allow chromosome names to start with a letter, so fix that, and sort by coordinates.
cat sample_ids_crown_rot.out_nodup.vcf | sed 's/LG//' | sed 's/Unknown/0/' | awk 'NR<3{print $0;next}{print $0| "sort -k1,2"}'  >sample_ids_crown_rot.out_nodup_fix.vcf
#Total genotyping rate is 0.481915.
#49992 variants removed due to minor allele threshold(s)
#45073 variants and 106 people pass filters and QC.
plink --vcf sample_ids_crown_rot.out_nodup_fix.vcf --maf 0.05 --recode --out sample_ids_crown_rot.out_nodup_fix_min05
#Change sex to male from Unknown
awk '{$5 = "1"; print}' sample_ids_crown_rot.out_nodup_fix_min05.ped > temp
mv temp sample_ids_crown_rot.out_nodup_fix_min05.ped 
#Substitute the missing phenotype values for mean crown rot scores.
python $scripts/add_phenotype_ped.py sample_ids_crown_rot.out_nodup_fix_min05.ped crown_rot_scores.txt >sample_ids_crown_rot.out_nodup_fix_min05_pheno.ped