#!/bin/bash

scripts=/home/sobczm/bin/popgen
#1st argument: Whole-genome GFF file.
gff=$1
filename=$(basename "$gff")
#First sort the gff file by chromosomes
sort -n -k 1 $gff >${filename%.gff*}_sorted.gff
#And then split
python $scripts/summary_stats/split_gff.py ${filename%.gff*}_sorted.gff
