#!/bin/bash
scripts=/home/sobczm/bin/popgen/snp
input=/home/sobczm/popgen/snp/snp_chip/crown_rot_gwas/plate/strawberry_db
#Various checks of the ananassa snp chip database - QC. 
cd $input

#Dump all the relevant strawberry_sample tables to files to be analyzed.
a="SELECT id, pipeline_id, sample_id, alias_id, clustering, genotype FROM genotype"
echo $a | mysql -u marias -h mongo -D strawberry_samples -p$(cat /home/sobczm/bin/mysql_sample_database/login) >genotype
b="SELECT id, platform, probe_id, snp_id, marker_id FROM alias"
echo $b | mysql -u marias -h mongo -D strawberry_samples -p$(cat /home/sobczm/bin/mysql_sample_database/login) >alias
c="SELECT id, name, ref, alt, upstream, downstream FROM marker"
echo $c | mysql -u marias -h mongo -D strawberry_samples -p$(cat /home/sobczm/bin/mysql_sample_database/login) >marker
d="SELECT marker_id, tag FROM marker_tag"
echo $d | mysql -u marias -h mongo -D strawberry_samples -p$(cat /home/sobczm/bin/mysql_sample_database/login) >marker_tag
e="SELECT id, clone_id, file, path, type, date, batch FROM sample"
echo $e | mysql -u marias -h mongo -D strawberry_samples -p$(cat /home/sobczm/bin/mysql_sample_database/login) >sample
#Print the genotype QC table for all samples.
#Take a given marker_id and print the number of genotypes available for a given marker_id on istraw35 and istraw90, and alias_id, Affx codes associated with it on each chip.
python $scripts/db_qc.py alias genotype sample >qc_table
Rscript --vanilla $scripts/analyze_db.R
#Print my list of multiform and secondaryform markers.
#multiform means "Affx name with more than one marker id associated", secondaryform mean "Affx name with more than one marker id associated, *this* form does not match any istraw35 marker"
python $scripts/db_qc2.py alias 
#Distribution of number of probes per the same Affx ID in istraw35 and istraw90.
python $scripts/db_qc3.py alias 
#Output:
## istraw35
#1    30014
#2    4246

## istraw90
#1    67561
#2    25736
#4    4
#6    347
#12    1414

#For cases, where istraw90 and istraw35 both contain one unique Affx, check if they match - if they don't print them out. Also, check the total number of matches.
python $scripts/db_qc4.py alias 

#For istraw35/90 multiforms separately, need to check that no duplicates linking the same marker Affx id to the same alleles with the same flanking sequence (but MAYBE use different probes) and the same marker id - highlight them, and choose one with more genotypes available to proceed in each chip. 
python $scripts/db_qc5.py alias marker
#Positive control test using fake data to check if script works correctly, as obtained no dupliactes
python $scripts/db_qc5.py alias_test marker_test

#Extract all possible samples, including duplicates in the spreadsheet.
python $scripts/ananassa_genotypes_db.py all_everbearer_samples.txt all_everbearer_samples.out
#All th samples in the entire db. 
python $scripts/ananassa_genotypes_db.py all_samples.txt all_samples.out

python $scripts/ananassa_genotypes_vcf.py all_everbearer_samples.out istraw90_vesca_v2.0_snp_positions.gff3
python $scripts/ananassa_genotypes_vcf.py all_samples.out istraw90_vesca_v2.0_snp_positions.gff3
#How many SNPs left if only 5% missing data allowed?
vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf all_everbearer_samples.out.vcf --max-missing 0.95 --recode --out all_everbearer_sample_max005
#kept 16419 out of a possible 95331 Sites
$vcftools/vcftools --vcf all_samples.out.vcf --max-missing 0.95 --recode --out all_samples_max005
#After filtering, kept 12932 out of a possible 95581 Sites

#Sort the resulting VCF files and fix so that compatible with Plink.
for a in *.vcf
do
$vcftools/vcf-sort $a | sed 's/LG//' | sed 's/Unknown/0/' >${a%.vcf}_sorted.vcf
done
#Calculate kinship in PLINK.
plink --vcf all_everbearer_sample_max005.recode_sorted.vcf --genome --out all_everbearer_kinship
cat all_everbearer_kinship.genome | awk '{$1=$1;print}' OFS='\t' >temp
mv temp all_everbearer_kinship.genome
plink --vcf all_samples_max005.recode_sorted.vcf --genome --out all_samples_kinship
cat all_samples_kinship.genome | awk '{$1=$1;print}' OFS='\t' >temp
mv temp all_samples_kinship.genome

Rscript --vanilla $scripts/plot_relationship.R 
Rscript --vanilla $scripts/plot_relationship2.R 

#Filter the rows where clone ids match for both samples being compared.
head all_samples_kinship_all.genome

#Plot my distance matrix
$scripts/similarity_percentage.py all_everbearer_sample_max005.recode_sorted.vcf
#Visualise the output as heatmap and clustering dendrogram
Rscript --vanilla $scripts/distance_matrix.R all_everbearer_sample_max005.recode_sorted_distance.log