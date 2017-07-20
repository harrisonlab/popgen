#!/bin/bash
#Note: fastStructure allows quick evaluation of population structure in DIPLOID organisms.
input=/home/sobczm/popgen/snp/snp_chip/fast_structure_test
scripts=/home/sobczm/bin/popgen/snp
pgdspid=/home/sobczm/bin/PGDSpider_2.1.0.3

cd $input
#Prepare STRUCTURE input (PGDSpider does not work when wrapped up in a bash script, grrr)
#haploid (for diploid change the conversion script to vcf_to_structure_diploid.spid)
#!!!! Need to change the path to file with population definitions !!!
input_file=$input/test.vcf
#In this file, chromosome numbers are preceded by LG - need to get rid of this prefix, as plink
#only accepts integers in the chromosome field.
sed -i 's/LG//' $input_file

#Convert from VCF to Plink's PED format
plink=/home/vicker/programs/plink-1.90beta/plink
$plink --vcf $input_file --const-fid 0 --out ${input_file%.vcf}

#Add '--prior logistic' when invoking structure inside the script below to change prior for more sensitive but exponentially slower structure detection.
s=1
f=5
for i in $(seq $s $f) #input range of K values tested
do
qsub $scripts/sub_fast_structure.sh ${input_file%.vcf} $i
done

#Choosing model complexity (K) among all the Ks tested
structure=/home/sobczm/bin/fastStructure
python $structure/chooseK.py --input=${input_file%.vcf}

#Visualize the expected admixture proportions inferred by fastStructure with Distruct plots using the mean of the variational posterior distribution over admixture proportions.

#First, generate sample labels:
cut -f2 ${input_file%.vcf}.fam >${input_file%.vcf}.lab

for i in $(seq $s $f)
do 
python $structure/distruct.py -K $i --input=${input_file%.vcf} --output=${input_file%.vcf}_${i}.svg --title K$i --popfile=${input_file%.vcf}.lab
done