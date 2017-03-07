#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 2
#$ -l virtual_free=5G
#$ -l h=blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

###Warning! Following the inspection of scree plot after the first run, 
#the R script then requires manual changing of the final_k parameter 
#which is the k value best explaining the population structure.

input=$1
ploidy=$2

scripts=/home/sobczm/bin/popgen/summary_stats
Rscript --vanilla $scripts/run_pcadapt.R $input $ploidy
