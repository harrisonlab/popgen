#!/bin/bash
#To change in each analysis:

input=/home/sobczm/popgen/input/mappings
reference=Fus2_canu_contigs_unmasked.fa

filename=$(basename "$reference")
output="${filename%.*}.dict"

#Prepare genome reference indexes required by GATK
java -jar /home/sobczm/bin/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R=$input/$reference O=$input/$output
samtools faidx $input/$reference

#Move to the directory where the output of SNP calling should be placed
cd /home/sobczm/popgen/snp/snp_calling/multithreaded

#Start SNP calling with GATK
#The submission script required need to be custom in each analysis, depending on what samples are being analysed.
scripts=/home/sobczm/bin/popgen/snp
qsub $scripts/sub_fus_SNP_calling_multithreaded.sh
