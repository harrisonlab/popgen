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

mkdir $input/first_attempt
cp apple*centroids.fasta 
cd $input/first_attempt

#Extract gene ids of centroids (separate Li and Velasco into 2 diff files at this point again)
for fasta in apple*.fasta
do
grep "^>" $fasta | sed 's/>//g' | grep "^evm.model*" >${fasta%.fasta}_li.lst
grep "^>" $fasta | sed 's/>//g' | grep "^MDP" >${fasta%.fasta}_velasco.lst
done 

genome_in=/home/sobczm/popgen/renseq/apple/genome
Li_gff=$genome_in/gene/EVM.out.gff
#Extract the corresponding gff intervals for exons
#Li
for a in *_li.lst
do
grep -f $a $Li_gff | awk '$3=="CDS"' >${a%.lst}.gff
done
#Velasco
Velasco_gff=$genome_in/Velasco/Malus_x_domestica.v3.0.a1_v1_gene_alignemnt.gff3
for a in *_velasco.lst
do
grep -f $a $Velasco_gff >${a%.lst}.gff
done


#Extract exon fasta sequences based on gff intervals. Use the masked version of Li genome and first hard-mask repeats with repeats from the DB provided with velasco genome. 
Li_assembly=$genome_in/assembly/malus_plantunas_100x.fasta.masked
for b in *li.gff
do
bedtools getfasta -fi $Li_assembly -bed $b -fo ${b%.gff}.fasta
done

Velasco_assembly=$genome_in/Velasco/Malus_x_domestica.v3.0.a1_contigs.fasta
Velasco_bed=$genome_in/Velasco/Malus_x_domestica.v3.0.a1_repeats.txt
bedtools maskfasta -fi $Velasco_assembly -bed $Velasco_bed -fo ${Velasco_assembly%.fasta}_masked.fasta
#Get an error: terminate called after throwing an instance of 'std::out_of_range'
#Needs to pre-filter the bed intervals to make sure they are present in the input file.

for b in *velasco.gff
do
bedtools getfasta -fi $Velasco_assembly -bed $b -fo ${b%.gff}.fasta
done
#The Velasco sequences have been found to contain ambigious characters denoting 
#polymorphisms. Substitute with random nucleotide matching the regex. 

#Join Li and Velasco input files again and design baits towards exon sequences at 2x coverage

for a in *li.fasta
do
cat $a ${a%li.fasta}velasco.fasta > ${a%li.fasta}recombined.fasta
done

for b in *recombined.fasta
do
n=2
python $scripts/create_baits.py --inp $b --coverage $n --out ${b%.fasta}_baits.fasta
done