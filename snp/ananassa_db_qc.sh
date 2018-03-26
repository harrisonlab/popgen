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
f="SELECT id, name, mat_id, pat_id FROM clone"
echo $f | mysql -u marias -h mongo -D strawberry_samples -p$(cat /home/sobczm/bin/mysql_sample_database/login) >clone
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

#Filter the rows where clone ids match for both samples being compared (except for NULLS).
cat all_samples_kinship_all.genome | awk '$15==$16 {print $0}' | awk '$15 !=  "NULL" && $16 != "NULL" {print $0}' | sort -n -k12 >identical_cultivar_threshold

#Identify any potential duplicate cultivar IDs, ignoring cases when two NULL cultivar IDs present
cat all_samples_kinship_all.genome | awk '$15!=$16 {print $0}' | awk '$12 > 0.966733 {print $0}' >potential_duplicates
#And now, focus just on double NULLs.
cat all_samples_kinship_all.genome | awk '$15== "NULL" && $16=="NULL"{print $0}' | awk '$12 > 0.966733 {print $0}' >potential_duplicates_null
#Plot my distance matrix
$scripts/similarity_percentage.py all_everbearer_sample_max005.recode_sorted.vcf
#Visualise the output as heatmap and clustering dendrogram
Rscript --vanilla $scripts/distance_matrix.R all_everbearer_sample_max005.recode_sorted_distance.log

#Find the samples id of a diverse panel of named cultivars based on their cultivar id.
cut -f1 list_different_cultivars_id.txt | tr '\n' ','
h="SELECT id FROM sample WHERE clone_id in (1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,3332,3334,3337,3338,3339,3340,3341,3342,3344,3345,3346,3347,3348,3349,3351,3352,3353,3357,3358,3359,3361,3362,3363,3394,3395,3396,3397,3398,3400,3401,3402,3403,3404,3405,3406,3407,3408,3409,3410,3411,3412,3413,3414,3415,3416,3417,3418,3419,3423,3590,3591,3592,3593,3596,3600,3606,3611,3616,3619,3622,3625,3628,3641,3644,3647,3654,3655,3656,3660,3661,3669,3670,3680,3681,3695,3701,3708,3744,3745,3755);"
echo $h | mysql -u marias -h mongo -D strawberry_samples -p$(cat /home/sobczm/bin/mysql_sample_database/login) >samples_differen_cultivars_id.lst

fn=samples_differen_cultivars_id
python $scripts/ananassa_genotypes_db.py ${fn}.lst samples_differen_cultivars_id.out
#All th samples in the entire db. 
python $scripts/ananassa_genotypes_vcf.py ${fn}.out istraw90_vesca_v2.0_snp_positions.gff3

#Sort the resulting VCF files and fix so that compatible with Plink.
$vcftools/vcf-sort ${fn}.out.vcf | sed 's/LG//' | sed 's/Unknown/0/' >${fn}.out_sorted.vcf
$vcftools/vcftools --vcf ${fn}.out_sorted.vcf --max-missing 0.95 --recode --out ${fn}_sorted_max_005

#Calculate kinship in PLINK.
plink --vcf ${fn}_sorted_max_005.recode.vcf --genome --out ${fn}_kinship
cat ${fn}_kinship.genome | awk '{$1=$1;print}' OFS='\t' >temp
mv temp ${fn}_kinship.genome 
Rscript --vanilla $scripts/plot_relationship3.R 