#!/bin/bash

input=/home/sobczm/popgen/renseq/input/transcriptomes
scripts=/home/sobczm/bin/popgen/renseq

#Annotate the four unmasked transcriptomes with NLR Parser:
#Run NLR Parser behind a screen as java config is messed up on SGE
cd $input/KIM/
sh $scripts/sub_nlrparser.sh GBRQ01_1_fsa_nt_combined_kim.fasta
cd $input/NZ/
sh $scripts/sub_nlrparser.sh GBGJ01_1_fsa_nt_nz.fasta
cd $input/RAJ/
sh $scripts/sub_nlrparser.sh GBJZ01_1_fsa_nt_raj.fasta
