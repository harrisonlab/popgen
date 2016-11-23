#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 5
#$ -l virtual_free=1G
#$ -l h=blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

#Look for R genes in a fasta file with protein sequences with RGAugury

input=$1
cpath=$PWD
temp_dir="$TMPDIR"
mkdir -p $temp_dir
cp $input $temp_dir
protein_fasta=$(basename "$input")
cd $temp_dir

rgaugury=/home/sobczm/bin/rgaugury/RGAugury.pl
perl $rgaugury -p $protein_fasta -c 5
cp -r * $cpath
rm -rf $temp_dir
