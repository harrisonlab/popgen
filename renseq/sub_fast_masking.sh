#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=1G
#$ -l h=blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

#Quick and reliable nucleotide low complexity sequence masking.
#For protein use seg algorithm.

fasta=$1

filename=$(basename "$fasta")
name=${filename%.*}

usearch=/home/sobczm/bin/usearch/usearch9.0.2132_i86linux32
$usearch -fastx_mask $fasta -qmask dust -fastaout "${name}"_masked.fasta -hardmask
