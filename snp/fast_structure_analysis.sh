#!/bin/bash
#Note: fastStructure allows quick evaluation of population structure in DIPLOID organisms.
input=/home/sobczm/popgen/snp/snp_chip/fast_structure_test
scripts=/home/sobczm/bin/popgen/snp

cd $input
input_file=$input/test2.vcf

#Convert from VCF to Plink's PED format
plink --allow-extra-chr --const-fid 0 --vcf $input_file --recode --make-bed --out ${input_file%.vcf} >${input_file%.vcf}.log

#Add '--prior logistic' when invoking structure inside the script below to change prior for more sensitive but exponentially slower structure detection. 
s=1 #minimum number of considered clusters
f=5 #maximum number of considered clusters
for i in $(seq $s $f) #input range of K values tested
do
qsub $scripts/sub_fast_structure.sh ${input_file%.vcf} $i
done

#Choosing model complexity (K) among all the Ks tested
structure=/home/sobczm/bin/fastStructure
python $structure/chooseK.py --input=${input_file%.vcf} >${input_file%.vcf}_K_choice

#Visualize the expected admixture proportions inferred by fastStructure with Distruct plots using the mean of the variational posterior distribution over admixture proportions.
#First, generate sample labels:
cut -f2 ${input_file%.vcf}.fam | cut -d" " -f2 >${input_file%.vcf}.lab

for i in $(seq $s $f)
do 
python $structure/distruct_mod.py -K $i --input=${input_file%.vcf} --output=${input_file%.vcf}_${i}.svg --title K$i --popfile=${input_file%.vcf}.lab
done