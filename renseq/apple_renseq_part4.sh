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