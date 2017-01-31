#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 5
#$ -l virtual_free=1G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace11.blacklace

#Look for R genes in a fasta file with protein sequences with RGAugury
#!!Notabene:
# I noticed a small bug in the pipeline in the script pfamscan.RGA.summary.pl line 54
#in the regular expression to match the TIR domain. As it can now match any domain name
#containing the "tir" element, it also incorrectly matches the Stirrup domain as a TIR domain.
#In this case, changing to case sensitive matching will solve the problem, I think,
#but not sure about the other domain names (LRR, #NB-ARC).

input=$1
rgaugury=/home/sobczm/bin/rgaugury/RGAugury.pl
perl $rgaugury -p $input -c 5
