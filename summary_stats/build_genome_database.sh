#!/bin/bash
#Builds a genome snpEff database used to annotate VCF files with SNP effects. Required input:
#1st argument: Unmasked genome FASTA file
#2nd argument: GFF file with gene annotations for the genome
#3rd argument: Unique genome name to be used to identify the genome when annotating VCF files

fasta=$1
gff=$2
genome_name=$3

snpeff=/home/sobczm/bin/snpEff

#Add the following lines to the section with databases:
echo "#" $genome_name genome >> $snpeff/snpEff.config
echo $genome_name.genome : $genome_name  >> $snpeff/snpEff.config

#Collect input files
#The folder name below should correspond to the genome name (Fus2v1.0.genome) specified above, without the ".genome" suffix.
#In the example case, Fus2v1.0.
mkdir $snpeff/data/$genome_name
cp $fasta $snpeff/data/$genome_name
cp $gff $snpeff/data/$genome_name

#Rename input files
gff_f=$(basename "$gff")
genome_f=$(basename "$fasta")

cd $snpeff/data/$genome_name
mv $gff_f genes.gff
mv $genome_f sequences.fa
#Build database using GFF3 annotation
java -jar $snpeff/snpEff.jar build -gff3 -v $genome_name
