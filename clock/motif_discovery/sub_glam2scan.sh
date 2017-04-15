#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

#Scan for a nucleotide motif with gaps (of non-constant length)
#Input:
#Fasta seqeuence to be scanned - 1st argument
#Glam motif definition file - 2nd argument
#How to prepare the Glam motif definiton file? Look at http://web.mit.edu/meme_v4.9.0/doc/glam2scan_man.html
#under "Motif format" heading

#Output: Directory named after input fasta file and file with the motif definition.

fasta=$1
motif=$2
ffilename=$(basename "$fasta")
mfilename=$(basename "$motif")
out=${mfilename%.*}_vs_${ffilename%.*}

meme=/home/sobczm/bin/meme_4.11.2/bin
$meme/glam2scan n -2 -n 100000 $motif $fasta -o $out
