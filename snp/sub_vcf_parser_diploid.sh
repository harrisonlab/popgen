#$ -S /bin/bash
#$ -cwd
#$ -pe smp 6
#$ -l h_vmem=4G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace

Filter the SNPs. Can allow a certain number of inds with N sites
vcf=$1

vcftools=/home/sobczm/bin/vcftools/bin
vcflib=/home/sobczm/bin/vcflib/bin

set -- "${1:-'default for 1'}" "${2:-'default 2'}" "${3:-'default 3'}"
