#!/bin/bash
scripts=/home/sobczm/bin/popgen/snp
input=/home/sobczm/popgen/snp/snp_chip/more_samples

#separate the analysis for samples from istraw90
#and istraw35. Include Holiday & Korona samples
mkdir istraw35 && cd istraw35
a="SELECT id FROM sample WHERE clone_id in (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,3332,3333,3334,3335,3336,3337,3338,3339,3340,3341,3342,3344,3345,3346,3347,3348,3349,3350,3351,3352,3353,3354,3355,3356,3357,3358,3359,3360,3361,3362,3363,3394,3395,3396,3397,3398,3400,3401,3402,3403,3404,3405,3406,3407,3408,3409,3410,3411,3412,3413,3414,3415,3416,3417,3418,3419,3423,3590,3591,3592,3593,3596,3600,3606,3611,3616,3619,3622,3625,3628,3641,3644,3647,3654,3655,3656,3657,3658,3659,3660,3661,3669,3670,3680,3681,3690,3695,3701,3708,3735,3744,3745,3749,3755,4134,4136,4158,4162,4186,4477,4479,4481,3721,3607,3612,3751,3731,3750,3739,3736,3646,3734,3722,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173) AND clone_id NOT IN (SELECT clone_id AS ci FROM clone_tag WHERE tag = 'rogue') AND type = 'istraw35';"

echo $a | mysql -u marias -h mongo -D strawberry_copy -p$(cat /home/sobczm/bin/mysql_sample_database/login)> samples_to_analyze.txt


mkdir istraw90 && cd istraw90
b="SELECT id FROM sample WHERE clone_id in (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,3332,3333,3334,3335,3336,3337,3338,3339,3340,3341,3342,3344,3345,3346,3347,3348,3349,3350,3351,3352,3353,3354,3355,3356,3357,3358,3359,3360,3361,3362,3363,3394,3395,3396,3397,3398,3400,3401,3402,3403,3404,3405,3406,3407,3408,3409,3410,3411,3412,3413,3414,3415,3416,3417,3418,3419,3423,3590,3591,3592,3593,3596,3600,3606,3611,3616,3619,3622,3625,3628,3641,3644,3647,3654,3655,3656,3657,3658,3659,3660,3661,3669,3670,3680,3681,3690,3695,3701,3708,3735,3744,3745,3749,3755,4134,4136,4158,4162,4186,4477,4479,4481,3721,3607,3612,3751,3731,3750,3739,3736,3646,3734,3722,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173) AND clone_id NOT IN (SELECT clone_id AS ci FROM clone_tag WHERE tag = 'rogue') AND type = 'istraw90';"
echo $b | mysql -u marias -h mongo -D strawberry_copy -p$(cat /home/sobczm/bin/mysql_sample_database/login)> samples_to_analyze.txt

#Execute the below in each of the folders above.
python $scripts/ananassa_genotypes_db.py samples_to_analyze.txt samples_to_analyze.out

#Output the genotypes in the VCF format with locations substituted  according to map positions relative to vesca 1.1. 
python $scripts/ananassa_genotypes_vcf.py samples_to_analyze.out istraw90_vesca_v1.1_snp_positions.gff3

c="SELECT DISTINCT(id), clone_id FROM sample;"
echo $c | mysql -u marias -h mongo -D strawberry_copy -p$(cat /home/sobczm/bin/mysql_sample_database/login)> sample_clone_id.txt

python $scripts/substitute_sample_names.py samples_to_analyze.out.vcf cultivar_names.txt sample_clone_id.txt

#Sort the VCF file.
cat samples_to_analyze.out_new_names.vcf | perl /home/sobczm/bin/vcftools/bin/vcf-sort > samples_to_analyze.out_new_names_sorted.vcf

#Remove any SNPs with missing data
vcftools=/home/sobczm/bin/vcftools/bin
filename=samples_to_analyze.out_new_names_sorted.vcf 
$vcftools/vcftools --vcf $filename --max-missing 1 --recode --out ${filename%.vcf}_nomissing

#Carry out PCA and plot the results
#First just get only sample num IDs.
sed 's/_\w\+//g' samples_to_analyze.out_new_names_sorted_nomissing.recode.vcf >samples_to_analyze.out_new_names_sorted_nomissing.recode.num.vcf
Rscript --vanilla $scripts/pca.R samples_to_analyze.out_new_names_sorted_nomissing.recode.num.vcf
#Eliminate Dover 474 Sample (outlier)

#Make an NJ tree 
$scripts/nj_tree.sh samples_to_analyze.out_new_names_sorted_nomissing.recode.vcf 2

#Calculate the index for percentage of shared SNP alleles between the individs.
$scripts/similarity_percentage.py samples_to_analyze.out_new_names_sorted_nomissing.recode.vcf
#Visualise the output as heatmap and clustering dendrogram
Rscript --vanilla $scripts/distance_matrix.R samples_to_analyze.out_new_names_sorted_nomissing.recode_distance.log

#Remove the outlier Dover 474 sample
input_file=samples_to_analyze.out_new_names_sorted.vcf
cut -d$'\t' -f-204,206- $input_file >samples_to_analyze.out_new_names_sorted_no474.vcf

#And carry out the PCA, NJ, distance matrix
vcftools=/home/sobczm/bin/vcftools/bin
filename=samples_to_analyze.out_new_names_sorted_no474.vcf
$vcftools/vcftools --vcf $filename --max-missing 1 --recode --out ${filename%.vcf}_nomissing

sed 's/_\w\+//g' samples_to_analyze.out_new_names_sorted_no474_nomissing.recode.vcf >samples_to_analyze.out_new_names_sorted_no474_nomissing.recode.num.vcf
Rscript --vanilla $scripts/pca.R samples_to_analyze.out_new_names_sorted_no474_nomissing.recode.num.vcf

#Make an NJ tree 
$scripts/nj_tree.sh samples_to_analyze.out_new_names_sorted_no474_nomissing.recode.vcf 2

#Calculate the index for percentage of shared SNP alleles between the individs.
$scripts/similarity_percentage.py samples_to_analyze.out_new_names_sorted_no474_nomissing.recode.vcf 
#Visualise the output as heatmap and clustering dendrogram
Rscript --vanilla $scripts/distance_matrix.R samples_to_analyze.out_new_names_sorted_no474_nomissing.recode_distance.log

#Calculate nucleotide diversity for different subsets.
vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf samples_to_analyze.out_new_names_sorted_filtered.recode.vcf --site-pi  --out nucleotide_diversity.all

#Analyze the origin of the Dover sample.
mkdir Dover && cd Dover
c="SELECT id FROM sample WHERE clone_id in (9,16,17,3229,3230,3231,3718,3721,3735) AND clone_id NOT IN (SELECT clone_id AS ci FROM clone_tag WHERE tag = 'rogue') AND type = 'istraw90';"
echo $c | mysql -u marias -h mongo -D strawberry_copy -p$(cat /home/sobczm/bin/mysql_sample_database/login)> samples_to_analyze.txt

