#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=0.1G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

#Input:
#Fasta file to be subsampled
#Number of random sequences to be printed
#example usage: qsub sub_fasta_subsample.sh test.fasta 100

fasta=$1
number=$2

ffilename=$(basename "$fasta")

RANDOM=`date +%N|sed s/...$//`

meme=/home/sobczm/bin/meme_4.11.2/bin/

$meme/fasta-subsample -seed $RANDOM $fasta $number >${ffilename%.*}_random_$number.fasta
