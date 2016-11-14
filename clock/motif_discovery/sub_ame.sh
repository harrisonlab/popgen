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
motif=$2
control=$3
ffilename=$(basename "$fasta")
mfilename=$(basename "$motif")
out=${mfilename%.*}_vs_${ffilename%.*}


meme=/home/sobczm/bin/meme_4.11.2/bin
$meme/ame --oc $out --control $control --pvalue-threshold 0.05 $fasta $motif 
