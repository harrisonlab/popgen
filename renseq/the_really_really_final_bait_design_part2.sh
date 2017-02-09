#!/bin/bash
input=/home/sobczm/popgen/renseq/input/transcriptomes
scripts=/home/sobczm/bin/popgen/renseq

#Separate the gene list into those on the positive and negative strands.

names=( "cornell" "han" "liu" "sun" "h6" "sp3b" "brian" "maria" "raj" "kim" "nz" )
for name in "${names[@]}"
do
cd $input/rgaugury/$name/lists
for a in *.strand
do
awk -F $"\t" '$2=="+" {print $1}' $a >"${a%.*}_pos.txt"
awk -F $"\t" '$2=="-" {print $1}' $a >"${a%.*}_neg.txt"
done
done

#Extract select contigs in each assembly
#genes on positive strand
for n in "${names[@]}"
do
cd $input/rgaugury/$n/lists
#Plus strand contigs
for list in *_pos.txt
do
python $scripts/keep_list_genes2.py $list $input/ultimate_baits/assemblies/${n}_assembly.fasta No
done
#Negative strand contigs - reverse complement
for list in *_neg.txt
do
python $scripts/keep_list_genes2.py $list $input/ultimate_baits/assemblies/${n}_assembly.fasta Yes
done
done

#Concatentate files with neg and pos strand contigs
for n in "${names[@]}"
do
cd $input/rgaugury/$n/lists
for fasta in *_pos.fasta
do
cat $fasta >>$input/really_really_final_baits/${fasta%_pos.fasta}.fasta
cat ${fasta%_pos.fasta}_neg.fasta >>$input/really_really_final_baits/${fasta%_pos.fasta}.fasta
done
done