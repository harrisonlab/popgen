#$ -S /bin/bash
#$ -cwd
#$ -pe smp 4
#$ -l h_vmem=1G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

#Fimo analysis to scan the sequences for a single motif.
#Input: 1st argument - a MEME Motif file (produced by iupac2meme)
#2nd argument - FASTA sequences to be scanned
#3rd argument - output dir

motif=$1
fasta=$2
out=$3

meme=/home/sobczm/bin/meme_4.11.2/bin
$meme/fimo --oc $out $motif $fasta
