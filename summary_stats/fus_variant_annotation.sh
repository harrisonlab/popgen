#!/bin/bash
#Copy input for the analysis into a new directory
input=/home/sobczm/popgen/summary_stats
snpeff=/home/sobczm/bin/snpEff
scripts=/home/sobczm/bin/popgen

#All individuals
cp /home/sobczm/popgen/snp/snp_calling/multithreaded/Fus2_canu_contigs_unmasked.vcf $input
cp /home/sobczm/popgen/snp/snp_calling/multithreaded/Fus2_canu_contigs_unmasked_filtered.vcf $input
cp /home/sobczm/popgen/input/Fus2_canu_new/final/Fus2_canu_contigs_unmasked.fa $input
cp /home/sobczm/popgen/input/Fus2_canu_new/final/Fus2_final_genes_appended.gff3 $input
cd $input
#Create additional subsets of VCF files with reduced number of individuals
vcftools=/home/sobczm/bin/vcftools/bin
#All, Without A13
$vcftools/vcftools --remove-indv FOCA13 --recode --recode-INFO-all --vcf Fus2_canu_contigs_unmasked.vcf --out noA13
$scripts/snp/vcf_parser_haploid.py --i noA13.recode.vcf
#Only pathogens, without A13
$vcftools/vcftools --remove-indv FOCA13 --remove-indv FOCA1-2 --remove-indv FOCD2 --remove-indv FOCA28 --remove-indv FOCCB3 --remove-indv FOCHB6 --remove-indv FOCPG --recode --recode-INFO-all --vcf Fus2_canu_contigs_unmasked.vcf --out patho
$scripts/snp/vcf_parser_haploid.py --i patho.recode.vcf
#Only non-pathogens, without A13
$vcftools/vcftools --remove-indv FOCA13 --remove-indv FOC55 --remove-indv FOCFus2 --remove-indv FOCA23 --remove-indv FOC125 --recode --recode-INFO-all --vcf Fus2_canu_contigs_unmasked.vcf --out non-patho
$scripts/snp/vcf_parser_haploid.py --i non-patho.recode.vcf

##Create custom SnpEff genome database
nano $snpeff/snpEff.config
#Add the following lines to the section with databases:
```
#---
# EMR Databases
#----
# Fus2 genome
Fus2v1.0.genome : Fus2
```
#Collect input files
mkdir $snpeff/data/Fus2v1.0
cp Fus2_canu_contigs_unmasked.fa $snpeff/data/Fus2v1.0
cp Fus2_final_genes_appended.gff3 $snpeff/data/Fus2v1.0
#Rename input files
cd $snpeff/data/Fus2v1.0
mv Fus2_final_genes_appended.gff3 genes.gff
mv Fus2_canu_contigs_unmasked.fa sequences.fa
#Build database using GFF3 annotation
java -jar $snpeff/snpEff.jar build -gff3 -v Fus2v1.0


#Annotate VCF files
cd $input
for a in *filtered.vcf
do
filename=$(basename "$a")
java -Xmx4g -jar $snpeff/snpEff.jar -v -ud 0 Fus2v1.0 $a > ${filename%.vcf}_annotated.vcf
mv snpEff_genes.txt snpEff_genes_${filename%.vcf}.txt
mv snpEff_summary.html  snpEff_summary__${filename%.vcf}.html
done
#Filter VCF files
for a in *annotated.vcf
do
filename=$(basename "$a")
#genic
java -jar $snpeff/SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') || (ANN[*].EFFECT has 'nonsense_variant') || (ANN[*].EFFECT has 'synonymous_variant') || (ANN[*].EFFECT has 'intron_variant') || (ANN[*].EFFECT has '5_prime_UTR_variant') || (ANN[*].EFFECT has '3_prime_UTR_variant')" $a > ${filename%.vcf}_gene.vcf
#coding
java -jar $snpeff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'synonymous_variant')" $a > ${filename%.vcf}_coding.vcf
#non-synonymous
java -jar $snpeff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant')" $a > ${filename%.vcf}_nonsyn.vcf
#synonymous
java -jar $snpeff/SnpSift.jar filter "(ANN[0].EFFECT has 'synonymous_variant')" $a > ${filename%.vcf}_syn.vcf
#silent
python $scripts/summary_stats/parse_snpeff_synonymous.py ${filename%.vcf}_syn.vcf
done
