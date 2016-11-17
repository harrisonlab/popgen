#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=0.1G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

#AME analysis to test if the select sequences are enriched for a specific motif
#Input: 1st argument - FASTA sequences to be scanned
#2nd argument - a MEME Motif file (produced by iupac2meme)
#3rd argument - FASTA sequences used as random control

fasta=$1
control=$2
ffilename=$(basename "$fasta")
out=${ffilename%.*}_glam

meme=/home/sobczm/bin/meme_4.11.2/bin
$meme/glam2 n -O $out -2 -n 1000 $fasta
#$meme/glam2 n -O $out -2 -n 1000000 $fasta
