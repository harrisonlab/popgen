#!/bin/bash

#Takes in a set of genome-wide CDS sequences in a FASTA file, filters them to remove those
#with premature stop codons, without a proper start codon at the beginning and stop codon at the end

#Subsequently, codon usage analysis is carried out with codonw.

f=$1

scripts=/home/sobczm/bin/popgen/codon
cw=/home/sobczm/popgen/codon/codonW

$scripts/check_proper_cds.py $f

filename=${f%.*}_pass.fasta

mkdir ${f%.*}
cp $filename ${f%.*}
cd ${f%.*}
# -cutot: tabulation of total codon usage
$cw/codonw $filename -nomenu -silent -cutot
mv ${filename%.*}.blk ${filename%.*}.cutot
# -cutab: tabulation of codon usage by gene
$cw/codonw $filename -nomenu -silent -cutab
mv ${filename%.*}.blk ${filename%.*}.cutab
# Correspondence analysis
# coa_num: percentage of genes taken from each tail of ENc distribution
$cw/codonw $filename -nomenu -silent -coa_cu -coa_num 5%
# Calculate all codon bias indices
# Use "cai.coa", "cbi.coa" and "fop.coa".  generated during correspondence
# Analysis to calculate the indices CAI, CBI and Fop
$cw/codonw $filename -all_indices -nomenu -silent -fop_file fop.coa -cai_file cai.coa -cbi_file cbi.coa
