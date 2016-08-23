#!/bin/bash

input=/home/sobczm/popgen/input/mappings
reference=Fus2_canu_contigs_unmasked.fa
filename=$(basename "$reference")
output="${filename%.*}.dict" 

#Prepare genome reference indexes
java -jar /home/sobczm/bin/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R=$input/$reference O=$input/$output
samtools faidx $input/$reference

cd /home/sobczm/popgen/snp/snp_calling
scripts=/home/sobczm/bin/scripts

qsub $scripts/sub_fus_SNP_calling.sh
