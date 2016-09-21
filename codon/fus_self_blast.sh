#!/bin/bash
input=/home/sobczm/popgen/input
wdir=/home/sobczm/popgen/codon/blast
scripts=/home/sobczm/bin/popgen/codon

#Copy input files
cd $wdir
cp $input/125/final/125_final_genes_combined.cdna.fasta $wdir
cp $input/55/final/55_final_genes_combined.cdna.fasta $wdir
cp $input/A1-2/final/A1-2_final_genes_combined.cdna.fasta $wdir
cp $input/A13/final/A13_final_genes_combined.cdna.fasta $wdir
cp $input/A23/final/A23_final_genes_combined.cdna.fasta $wdir
cp $input/A28/final/A28_final_genes_combined.cdna.fasta $wdir
cp $input/CB3/final/CB3_final_genes_combined.cdna.fasta $wdir
cp $input/D2/final/D2_final_genes_combined.cdna.fasta $wdir
cp $input/Fus2_canu_new/final/Fus2_final_genes_combined.cdna.fasta $wdir
cp $input/HB6/final/HB6_final_genes_combined.cdna.fasta $wdir
cp $input/PG/final/PG_final_genes_combined.cdna.fasta $wdir
cp $input/proliferatum/final/proliferatum_final_genes_combined.cdna.fasta $wdir

#Retain only one gene transcript per gene
for f in *.fasta
do
$scripts/keep_one_gene.py $f
done

#Create nucleotide databases
#Run assembly self-on-self BLAST
for assembly in *one.fasta
do
makeblastdb -in $assembly -input_type fasta -dbtype nucl \
-title "${assembly%.*}"_nucl.db -parse_seqids -out "${assembly%.*}"_nucl.db
qsub $scripts/sub_self_blast.sh $assembly "${assembly%.*}"_nucl.db
done
