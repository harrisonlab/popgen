#!/bin/bash
 
###Preparing alignments and finding best-fit nucleotide sequence evolution models

path=/home/sobczm/popgen
scripts=/home/sobczm/bin/scripts

##MAFFT (alignments)

cd $path/busco
cp *.fasta $path/busco_alignments 
cd $path/busco_alignments

for f in *.fasta; do qsub $scripts/sub_mafft_alignment.sh $f; done

## Identify genes with high nucleotide diversity

python $scripts/calculate_nucleotide_diversity.py
