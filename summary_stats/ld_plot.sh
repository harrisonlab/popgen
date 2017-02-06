#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace

input_ld=$1
scripts=/home/sobczm/bin/popgen/summary_stats

sed "s/\s*$//" $input_ld >${input_ld}_plotld
Rscript --vanilla $scripts/ld_plot.R ${input_ld}_plotld