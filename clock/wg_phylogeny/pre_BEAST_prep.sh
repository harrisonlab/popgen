#!/bin/bash

###Preparing alignments and finding best-fit nucleotide sequence evolution models

path=/home/sobczm/popgen/phylogenetics/clock/CDS_genomes
scripts=/home/sobczm/bin/popgen/phylogenetics

##MAFFT (make alignments)
cd $path
mkdir busco_alignments
mv *.fasta $path/busco_alignments
cd $path/busco_alignments

qsub $scripts/sub_mafft_alignment.sh $f

## Identify genes with high nucleotide diversity and average number of pairwise differences, medium number of segregating sites
## (avoid alignments with low homology and lots of phylogenetically uninformative singletons)

python $scripts/calculate_nucleotide_diversity.py "*aligned.fasta"

mkdir -d $path/beast_runs/results
mv sequence_stats.txt excel_stats.txt $path/beast_runs/results
