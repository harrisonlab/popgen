#!/bin/bash
input=/home/sobczm/popgen/clock/DNA_genomes/promoters/extended
scripts=/home/sobczm/bin/popgen/clock/motif_discovery
meme=/home/sobczm/bin/meme_4.11.2/bin
weeder=/home/sobczm/bin/weeder/weeder2

mkdir $input/ace
mkdir $input/cbox
################Motif discovery testing
# A) ACE motif containing (TCTTGGCA)
# B) Clockbox motif containing (CGAT(N)CCGCT)

#Fish out the homologs of A) and B) motif-containing Neurospora genes.
#To accomplish that, first parse the input spreadsheet with protein names
#to produce a custom list for each species
cp ace_all_genes.txt $input/ace
cp cbox_all_genes.txt $input/cbox
python $scripts/prepare_gene_list.py $input/ace/ace_all_genes.txt ace
for a in *ace; do sed -i '/^$/d' $a; done
python $scripts/prepare_gene_list.py $input/cbox/cbox_all_genes.txt cbox
for a in *cbox; do sed -i '/^$/d' $a; done

#for each genome extract promoter sequences for given genes in the list.
#ACE genes and Clockbox genes
sh $scripts/extract_ace_cbox_prom.sh

#First, use the degenerate A) and B) motifs identified in the previous analysis
#to test for enrichment amongst the putative circadian genes identified in Nc microarray.
#Second of all, carry out de novo motif discovery amongst the genes in each species' set

#Count the number of sequences in each file and create an appropriate background
#control file with the same number of random promoters.
for file in $input/ace/*filtered.fasta
do
lines=`wc -l $file | cut -f1 -d' '`
ngenes="$((lines/2))"
ffilename=$(basename "$file")
assembly=${ffilename%_filtered.fasta}.fasta
qsub $scripts/sub_fasta_subsample.sh $input/$assembly $ngenes
done

for file in $input/cbox/*filtered.fasta
do
lines=`wc -l $file | cut -f1 -d' '`
ngenes="$((lines/2))"
ffilename=$(basename "$file")
assembly=${ffilename%_filtered.fasta}.fasta
qsub $scripts/sub_fasta_subsample.sh $input/$assembly $ngenes
done

#DREME and GLAM2
for control in $input/cbox/*random*.fasta
do
file=${control%_random*}_filtered.fasta
qsub $scripts/sub_dreme.sh $file $control
qsub $scripts/sub_glam.sh $file
done

for control in $input/ace/*random*.fasta
do
file=${control%_random*}_filtered.fasta
qsub $scripts/sub_dreme.sh $file $control
qsub $scripts/sub_glam.sh $file
done
