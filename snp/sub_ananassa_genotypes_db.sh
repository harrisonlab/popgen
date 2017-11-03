#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=16G
#$ -l h=blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace
input=$1
output=$2

/home/sobczm/bin/popgen/snp/ananassa_genotypes_db.py $input $output