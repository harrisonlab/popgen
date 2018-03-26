#!/bin/bash
scripts=/home/sobczm/bin/popgen/renseq
input=/home/sobczm/popgen/renseq/apple_2017

#Extract the corresponding gff intervals for exons
#Li
for a in *no_transposon
do
grep -f $a gene_models_20170612.gff3 | awk '$3=="exon"' >${a%.lst}.gff
done

#Extract exon fasta sequences based on gff intervals. First hard-mask repeats with repeats from the DB provided with the genome. 
assembly=GDDH13_1-1_formatted.fasta
repeats=GDDH13_1-1_TE.gff3
bedtools maskfasta -fi $assembly -bed $repeats -fo ${assembly%.fasta}_masked.fasta

for b in *no_transposon.gff
do
bedtools getfasta -fi ${assembly%.fasta}_masked.fasta -bed $b -fo ${b%.gff}.fasta
done

#Design baits towards exon sequences at 2x coverage
for b in *no_transposon.fasta
do
n=2
python $scripts/create_baits.py --inp $b --coverage $n --out ${b%.fasta}_baits.fasta
done
