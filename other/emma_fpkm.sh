#!/bin/bash

input=/home/sobczm/popgen/other/emma_fpkm
#CDS lengths to use for FPKM calculations
python /home/sobczm/bin/popgen/renseq/write_seq_length.py Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all.fa
#Remove the .t1 etc. transcript suffix in the file with CDS lengths info
#There is only one transcript per gene, so no need for further manipulation
sed -i 's/.t[0-9]//' Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all_lengths.txt 
#Add headers to the length table and fix the header in the count table.
sed -i '1s/^/Gene\tLength\n/' Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all_lengths.txt
#Calculate FPKM with a new script
scripts=/home/sobczm/bin/popgen/rnaseq

python $scripts/calculate_fpkm.py countData_53 Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all_lengths.txt
 
python $scripts/calculate_fpkm.py countData_08 Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all_lengths.txt
 