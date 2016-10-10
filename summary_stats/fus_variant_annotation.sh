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
vcflib=/home/sobczm/bin/vcflib/bin

#All, without A13
$vcflib/vcfremovesamples Fus2_canu_contigs_unmasked.vcf FOCA13 >Fus2_canu_contigs_unmasked_noA13.vcf
#Remove monomorphic sites (minor allele count minimum 1)
$scripts/snp/vcf_parser_haploid.py --i Fus2_canu_contigs_unmasked_noA13.vcf
$vcftools/vcftools --vcf Fus2_canu_contigs_unmasked_noA13_filtered.vcf --mac 1 --recode --out Fus2_canu_contigs_unmasked_noA13_filtered

#Only pathogens, without A13
$vcflib/vcfremovesamples Fus2_canu_contigs_unmasked.vcf FOCA13 FOCA1-2 FOCD2 FOCA28 FOCCB3 FOCHB6 FOCPG >Fus2_canu_contigs_unmasked_patho.vcf
$scripts/snp/vcf_parser_haploid.py --i Fus2_canu_contigs_unmasked_patho.vcf
$vcftools/vcftools --vcf Fus2_canu_contigs_unmasked_patho_filtered.vcf --mac 1 --recode --out Fus2_canu_contigs_unmasked_patho_filtered

#Only pathogens, without A13 and 55
$vcflib/vcfremovesamples Fus2_canu_contigs_unmasked.vcf FOCA13 FOCA1-2 FOCD2 FOCA28 FOCCB3 FOCHB6 FOCPG FOC55 >Fus2_canu_contigs_unmasked_patho_no55.vcf
$scripts/snp/vcf_parser_haploid.py --i Fus2_canu_contigs_unmasked_patho_no55.vcf
$vcftools/vcftools --vcf Fus2_canu_contigs_unmasked_patho_no55_filtered.vcf --mac 1 --recode --out Fus2_canu_contigs_unmasked_patho_no55_filtered

#Only non-pathogens, without A13
$vcflib/vcfremovesamples Fus2_canu_contigs_unmasked.vcf FOCA13 FOC55 FOCFus2 FOCA23 FOC125 >Fus2_canu_contigs_unmasked_non-patho.vcf
$scripts/snp/vcf_parser_haploid.py --i Fus2_canu_contigs_unmasked_non-patho.vcf
$vcftools/vcftools --vcf Fus2_canu_contigs_unmasked_non-patho_filtered.vcf --mac 1 --recode --out Fus2_canu_contigs_unmasked_non-patho_filtered

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
for a in *recode.vcf
do
filename=$(basename "$a")
java -Xmx4g -jar $snpeff/snpEff.jar -v -ud 0 Fus2v1.0 $a > ${filename%.vcf}_annotated.vcf
mv snpEff_genes.txt snpEff_genes_${filename%.vcf}.txt
mv snpEff_summary.html  snpEff_summary__${filename%.vcf}.html
done
###Further filter VCF files by SNP class
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


###Create FASTA alignment files containing only select subsets of SNPs
### From now onwards, analysing the dataset without A13.
cd $input/noA13
ref_genome=/home/sobczm/popgen/summary_stats/Fus2_canu_contigs_unmasked.fa
##all
python $scripts/summary_stats/vcf_to_fasta.py Fus2_canu_contigs_unmasked_noA13_filtered.recode_annotated.vcf $ref_genome 1
#Moving each subset of FASTA files into a separate dir.
mkdir all
mv *.fasta ./all
##coding
python $scripts/summary_stats/vcf_to_fasta.py Fus2_canu_contigs_unmasked_noA13_filtered.recode_annotated_coding.vcf $ref_genome 1
mkdir coding
mv *.fasta ./coding
##silent(four-fold degenerate)
python $scripts/summary_stats/vcf_to_fasta.py Fus2_canu_contigs_unmasked_noA13_filtered.recode_annotated_syn_silent.vcf $ref_genome 1
mkdir silent
mv *.fasta ./silent


###Split the GFF file into one contig --> one GFF file.
cd $input
mkdir gff
#First sort the gff file by chromosomes
sort -n -k 1 Fus2_final_genes_appended.gff3 >Fus2_final_genes_appended_sorted.gff3
#And then split
python $scripts/summary_stats/split_gff.py Fus2_final_genes_appended_sorted.gff3
mv *.gff ./gff
