#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=0.1G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

#Fimo analysis to scan the sequences for a single motif.without gaps
#Input: 1st argument - FASTA sequences to be scanned
#2nd argument - give a name for the motif you are scanning for (no spaces allowed)
#3nd argument - the 5'->3' sequence of the motif

#Output: Directory named after input fasta file and motif name,
#Motif definition file in the MEME format.

fasta=$1
motif_name=$2
motif=$3

ffilename=$(basename "$fasta")
out=${mfilename%.*}_vs_${motif_name}
meme=/home/sobczm/bin/meme_4.11.2/bin

$meme/iupac2meme $motif >${motif_name}.txt
$meme/fimo --oc $out ${motif_name}.txt $fasta
