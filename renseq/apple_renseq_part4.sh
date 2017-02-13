#!/bin/bash
scripts=/home/sobczm/bin/popgen/renseq
input=/home/sobczm/popgen/renseq/apple/baits

#Contig clustering 
cd $input
usearch=/home/sobczm/bin/usearch/usearch9.0.2132_i86linux32
for fasta in apple*.fasta
do
id=0.9
$usearch -cluster_fast $fasta -id $id -sort length -centroids ${fasta%.fasta}_${id}_centroids.fasta
done 

#Design baits towards gene sequences (including introns)
#at 4x coverage
#!!Test!!: just CDS sequences.
#First, hard-mask repetitive sequences in each centroid file
for fasta in *centroid*.fasta
do
$usearch -fastx_mask $fasta -qmask dust -fastaout "${fasta%.*}"_masked.fasta -hardmask
done

#The Velasco sequences have been found to contain ambigious characters
#Substitute at random


#Hard-mask repeats with repeats from the DB