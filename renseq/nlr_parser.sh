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

#Extract target putative R contigs ID. Exclude sequences with domain matches
#in plus and minus strands as likely result of merging of different transcripts
#during the assembly process.


#Extract the target sequences and reverse complement sequences on the negative strand.


#Design baits to the select transcripts
