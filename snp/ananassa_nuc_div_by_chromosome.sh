#!/bin/bash
scripts=/home/sobczm/bin/popgen/snp
input=/home/sobczm/popgen/snp/snp_chip/more_samples/detailed_nuc_diversity

#Extract samples from istraw90. 
b="SELECT id FROM sample WHERE clone_id in (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,3332,3333,3334,3335,3336,3337,3338,3339,3340,3341,3342,3344,3345,3346,3347,3348,3349,3350,3351,3352,3353,3354,3355,3356,3357,3358,3359,3360,3361,3362,3363,3394,3395,3396,3397,3398,3400,3401,3402,3403,3404,3405,3406,3407,3408,3409,3410,3411,3412,3413,3414,3415,3416,3417,3418,3419,3423,3590,3591,3592,3593,3596,3600,3606,3611,3616,3619,3622,3625,3628,3641,3644,3647,3654,3655,3656,3657,3658,3659,3660,3661,3669,3670,3680,3681,3690,3695,3701,3708,3735,3744,3745,3749,3755,4134,4136,4158,4162,4186,4477,4479,4481,3721,3607,3612,3751,3731,3750,3739,3736,3646,3734,3722) AND clone_id NOT IN (SELECT clone_id AS ci FROM clone_tag WHERE tag = 'rogue') AND type = 'istraw90';"
echo $b | mysql -u marias -h mongo -D strawberry_copy -p$(cat /home/sobczm/bin/mysql_sample_database/login)> samples_to_analyze.txt

#Extract genotypes for selected samples for downstream processing.
python $scripts/ananassa_genotypes_db.py samples_to_analyze.txt samples_to_analyze.out

#Output the genotypes in the VCF format with locations substituted 
#according to map positions relative to vesca 1.1. 
python $scripts/ananassa_genotypes_vcf.py samples_to_analyze.out istraw90_vesca_v1.1_snp_positions.gff3

#Split the VCF file by chromomsome.
python $scripts/split_vcf_by_chromosome.py samples_to_analyze.out.vcf vesca2consensus_map_noduplicates_2017-05-17.csv

#Calculate the nucleotide diversity by chromomsome.
vcftools=/home/sobczm/bin/vcftools/bin
for vcf in *.vcf
do
$vcftools/vcftools --vcf $vcf --site-pi --out ${vcf%.vcf}.pi
done

for pi in *.pi
do
Rscript --vanilla $scripts/summary_distribution.R $pi
done

#Collect the results in one file.
echo "Min.  1st_Qu. Median  Mean    3rd_Qu. Max." >pi_all_chroms.txt
for file in *.stat
do
echo $file >>pi_all_chroms.txt
cat $file >>pi_all_chroms.txt
done

#Convert spaces to tabs
cat pi_all_chroms.txt | awk '{$1=$1;print}' OFS='\t' 