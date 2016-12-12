#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=0.1G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

#GLAM analysis to try to identify GAPPED short motifs for which a given set of sequences is significantly enriched for
#Input: input FASTA sequences to be used for motif discovery

#Output:
#Folder (suffix "glam") containing sequence logos, text files with sequence motifs

fasta=$1
ffilename=$(basename "$fasta")
out=${ffilename%.*}_glam

meme=/home/sobczm/bin/meme_4.11.2/bin
$meme/glam2 n -O $out -2 -n 10000 $fasta
