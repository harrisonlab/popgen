#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=0.2G
#$ -l h=blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

#Look for R genes in a fasta file with protein sequences with RGAugury
#!!Notabene:
# I noticed a small bug in the pipeline in the script pfamscan.RGA.summary.pl line 54
#in the regular expression to match the TIR domain. As it can now match any domain name
#containing the "tir" element, it also incorrectly matches the Stirrup domain as a TIR domain.
#In this case, changing to case sensitive matching will solve the problem, I think,
#but not sure about the other domain names (LRR, NB-ARC).

input=$1
cpath=$PWD
temp_dir="$TMPDIR"
mkdir -p $temp_dir
cp $input $temp_dir
protein_fasta=$(basename "$input")
cd $temp_dir

rgaugury=/home/sobczm/bin/rgaugury2/RGAugury.pl
perl $rgaugury -p $protein_fasta -c 1
cp -r * $cpath
rm -rf $temp_dir
