#!/bin/bash
input=/home/sobczm/popgen/renseq/input/transcriptomes/plant_rgenes/common_onion/domains
scripts=/home/sobczm/bin/popgen/renseq

#Detect contigs containing two other domains in the PFAM output from plant_rgenes.sh
#for each transcriptome.
#PF05627	Cleavage site for pathogenic type III effector avirulence factor Avr
#PF05659    RPW8

cd $input
names=( "kim" "nz" "maria" "raj" "brian" "cornell" "h6" "han" "sp3b" "liu" "sun" )
for name in "${names[@]}"
do
cd $input/$name
grep PF05659 onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.out | awk -v \
OFS="\t" '$1=$1' >${name}_rpw8.txt
grep PF05627 onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.out | awk -v \
OFS="\t" '$1=$1' >${name}_Avr_cleaveage.txt
grep "NB-ARC" onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.out | sort | awk -v \
OFS="\t" '$1=$1' >${name}_NBS.txt
python $scripts/add_strand_column.py ${name}_rpw8.txt
python $scripts/add_strand_column.py ${name}_Avr_cleaveage.txt
python $scripts/add_strand_column.py ${name}_NBS.txt
done

#Select only good matches, with p-value < e-05
