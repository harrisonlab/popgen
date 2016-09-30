#!/bin/bash

wdir=/home/sobczm/popgen/codon/blast
input=/home/sobczm/popgen/input
scripts=/home/sobczm/bin/popgen/codon

#Copy over relevant GFF files
cp $input/125/final/*appended.gff3 $wdir
cp $input/55/final/*appended.gff3 $wdir
cp $input/A1-2/final/*appended.gff3 $wdir
cp $input/A13/final/*appended.gff3 $wdir
cp $input/A23/final/*appended.gff3 $wdir
cp $input/A28/final/*appended.gff3 $wdir
cp $input/CB3/final/*appended.gff3 $wdir
cp $input/D2/final/*appended.gff3 $wdir
cp $input/Fus2_canu_new/final/*appended.gff3 $wdir
cp $input/HB6/final/*appended.gff3 $wdir
cp $input/PG/final/*appended.gff3 $wdir
cp $input/proliferatum/final/*appended.gff3 $wdir


## Parse the GFF files and prepare annotation tables
for a in *_final_genes_appended.gff3
do
$scripts/cds_to_chromosome_coords.py $a ${a%%final_genes_appended.gff3} #trim the match from the end and extract prefix
done

## Filter BLAST matches by subject coverage and sequence Identity
for b in *.db
do
$scripts/filter_blast.py --i $b
done
