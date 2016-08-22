#!/bin/bash

input=/home/sobczm/popgen/input/mappings
scripts=/home/sobczm/bin/scripts

### Rename input files in each folder by prefixing with the strain ID
cd $input/125/
for filename in *; do mv "$filename" "125_$filename"; done;
cd $input/55/
for filename in *; do mv "$filename" "55_$filename"; done;
cd $input/A1-2/
for filename in *; do mv "$filename" "A1-2_$filename"; done;
cd $input/A13/
for filename in *; do mv "$filename" "A13_$filename"; done;
cd $input/A23/
for filename in *; do mv "$filename" "A23_$filename"; done;
cd $input/A28/
for filename in *; do mv "$filename" "A28_$filename"; done;
cd $input/CB3/
for filename in *; do mv "$filename" "CB3_$filename"; done;
cd $input/D2/
for filename in *; do mv "$filename" "D2_$filename"; done;
cd $input/Fus2/
for filename in *; do mv "$filename" "Fus2_$filename"; done;
cd $input/HB6/
for filename in *; do mv "$filename" "HB6_$filename"; done;
cd $input/PG/
for filename in *; do mv "$filename" "PG_$filename"; done;

### Remove multimapping reads, discordant reads. PCR and optical duplicates, and 
### add read group and sample name to each mapped read

qsub $scripts/sub_pre_snp_calling.sh $input/125/125_Fus2_canu_contigs_unmasked.fa_aligned.sam FOC125
qsub $scripts/sub_pre_snp_calling.sh $input/55/55_Fus2_canu_contigs_unmasked.fa_aligned.sam FOC55
qsub $scripts/sub_pre_snp_calling.sh $input/A1-2/A1-2_Fus2_canu_contigs_unmasked.fa_aligned.sam FOCA1-2
qsub $scripts/sub_pre_snp_calling.sh $input/A13/A13_Fus2_canu_contigs_unmasked.fa_aligned.sam FOCA13
qsub $scripts/sub_pre_snp_calling.sh $input/A23/A23_Fus2_canu_contigs_unmasked.fa_aligned.sam FOCA23
qsub $scripts/sub_pre_snp_calling.sh $input/A28/A28_Fus2_canu_contigs_unmasked.fa_aligned.sam FOCA28
qsub $scripts/sub_pre_snp_calling.sh $input/CB3/CB3_Fus2_canu_contigs_unmasked.fa_aligned.sam FOCCB3
qsub $scripts/sub_pre_snp_calling.sh $input/D2/D2_Fus2_canu_contigs_unmasked.fa_aligned.sam FOCD2
qsub $scripts/sub_pre_snp_calling.sh $input/Fus2/Fus2_Fus2_canu_contigs_unmasked.fa_aligned.sam FOCFus2
qsub $scripts/sub_pre_snp_calling.sh $input/HB6/HB6_Fus2_canu_contigs_unmasked.fa_aligned.sam FOCHB6
qsub $scripts/sub_pre_snp_calling.sh $input/PG/PG_Fus2_canu_contigs_unmasked.fa_aligned.sam FOCPG
