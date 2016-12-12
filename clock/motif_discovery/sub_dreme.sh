#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=0.1G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

#DREME analysis to try to identify UNGAPPED short motifs for which a given set of sequences is significantly enriched for
#Input:
#1st argument - input FASTA sequences to be used for motif discovery
#2nd argument - FASTA sequences to be used as background control (ie. not expected to
# contain a shared motif). They can be subsampled from e.g. file containing equivalent
#genomic sequences using sub_fasta_subsample.sh script.
#Output:
#Folder (suffix "dreme") containing sequence logos, text files with sequence motifs


fasta=$1
control=$2
ffilename=$(basename "$fasta")
out=${ffilename%.*}_dreme

meme=/home/sobczm/bin/meme_4.11.2/bin
$meme/dreme -oc $out -p $fasta -n $control -g 10000
