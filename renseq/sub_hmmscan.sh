#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 2
#$ -l virtual_free=1G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace11.blacklace

input_fasta=$1
filename=$(basename "$input_fasta")
name=${filename%.*}

pfam_scan.pl -e_seq 1 -e_dom 1 -as -outfile "${name}".out -cpu 2 \
-fasta $input_fasta -dir /home/sobczm/bin/hmmer-3.1b2/pfam30
