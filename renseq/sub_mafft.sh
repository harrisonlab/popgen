#$ -S /bin/bash
#$ -cwd
#$ -pe smp 6
#$ -l h_vmem=1G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace

input=$1
mafft_path=/home/sobczm/bin/mafft-7.222/bin/mafft

filename=$(basename "$input")
output="${filename%.*}_aligned.fasta"
$mafft_path --localpair --maxiterate 1000 $input >$output