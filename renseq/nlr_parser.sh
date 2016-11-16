#!/bin/bash
input=/home/sobczm/popgen/renseq/input/transcriptomes
scripts=/home/sobczm/bin/popgen/renseq

#Annotate the unmasked transcriptomes with NLR Parser:
#Run NLR Parser behind a screen as java config is messed up on SGE and sort output
cd $input/KIM/
sh $scripts/sub_nlrparser.sh GBRQ01_1_fsa_nt_combined_kim.fasta
sort -k 1 GBRQ01_1_fsa_nt_combined_kim_nlr.tsv >GBRQ01_1_fsa_nt_combined_kim_nlr_sorted.tsv
cd $input/NZ/
sh $scripts/sub_nlrparser.sh GBGJ01_1_fsa_nt_nz.fasta
sort -k 1 GBGJ01_1_fsa_nt_nz_nlr.tsv >GBGJ01_1_fsa_nt_nz_nlr_sorted.tsv
cd $input/RAJ/
sh $scripts/sub_nlrparser.sh GBJZ01_1_fsa_nt_raj.fasta
sort -k 1 GBJZ01_1_fsa_nt_raj_nlr.tsv >GBJZ01_1_fsa_nt_raj_nlr_sorted.tsv
cd $input/MARIA/
sh $scripts/sub_nlrparser.sh Trinity.fasta
sort -k 1 Trinity_nlr.tsv >Trinity_nlr_sorted.tsv
cd $input/BRIAN/
sh $scripts/sub_nlrparser.sh Trinity_11Oct2016.fasta
sort -k 1 Trinity_11Oct2016_nlr.tsv >Trinity_11Oct2016_nlr_sorted.tsv
