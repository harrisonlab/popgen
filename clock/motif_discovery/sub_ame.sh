#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=0.1G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

#AME analysis to test if the select sequences are enriched for a specific motif with fixed length
#Input:
#1st argument - a set of FASTA sequences to be scanned.
#2nd argument: a set of background sequences to be used for controls. Can be
#subsampled from e.g. file containing equivalent genomic sequences using
#sub_fasta_subsample.sh script.
#3rd argument - give a name for the motif you are scanning for (no spaces allowed)
#4th argument - the 5'->3' sequence of the motif

#Output: Directory named after input fasta file and motif name,
#Motif definition file in the MEME format.

fasta=$1
control=$2
motif_name=$3
motif=$4

ffilename=$(basename "$fasta")
out=${mfilename%.*}_vs_$motif_name
meme=/home/sobczm/bin/meme_4.11.2/bin

$meme/iupac2meme $motif >${motif_name}.txt
$meme/ame --oc $out --control $control --pvalue-threshold 0.05 $fasta ${motif_name}.txt
