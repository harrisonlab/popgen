#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=16G
#$ -l h=blacklace11.blacklace
input=$1
output=$2

/home/sobczm/bin/popgen/snp/ananassa_genotypes_db.py $input $output