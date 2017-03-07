#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 20
#$ -l virtual_free=1G
#$ -l h=blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

input=$1
cpath=$PWD

temp_dir="$TMPDIR"
mkdir -p $temp_dir

cp -r $input $temp_dir
input_f=$(basename "$input")
cd $temp_dir

bayescan=/home/sobczm/bin/bayescan2.1/binaries/bayescan_2.1
$bayescan -threads 16 -od ./ $input_f

rm $input_f
cp -r * $cpath
rm -rf $temp_dir

scripts=/home/sobczm/bin/popgen/summary_stats
cd $cpath
Rscript --vanilla $scripts/plot_bayescan.R ${input%vcf}_fst.txt
