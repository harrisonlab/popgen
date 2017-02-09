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

#Concatenate the contig sequence files for each gene family
cd $input/really_really_final_baits
for a in *.NBS.all.fasta
do
cat $a >>onion.NBS.all.fasta
done

for a in *.NBS.CN.fasta
do
cat $a >>onion.NBS.CN.fasta
done

for a in *.NBS.CNL.fasta
do
cat $a >>onion.NBS.CNL.fasta
done

for a in *.NBS.NBS.fasta
do
cat $a >>onion.NBS.NBS.fasta
done

for a in *.NBS.NL.fasta
do
cat $a >>onion.NBS.NL.fasta
done

for a in *.NBS.OTHER.fasta
do
cat $a >>onion.NBS.OTHER.fasta
done

for a in *.NBS.TN.fasta
do
cat $a >>onion.NBS.TN.fasta
done

for a in *.NBS.TNL.fasta
do
cat $a >>onion.NBS.TNL.fasta
done

for a in *.NBS.TX.fasta
do
cat $a >>onion.NBS.TX.fasta
done

for a in *.RLK.all.fasta
do
cat $a >>onion.RLK.all.fasta
done

for a in *.RLK.LRR.fasta
do
cat $a >>onion.RLK.LRR.fasta
done

for a in *.RLK.LysM.fasta
do
cat $a >>onion.RLK.LysM.fasta
done

for a in *.RLK.other.fasta
do
cat $a >>onion.RLK.other.fasta
done

for a in *.RLP.all.fasta
do
cat $a >>onion.RLP.all.fasta
done

for a in *.RLP.LRR.fasta
do
cat $a >>onion.RLP.LRR.fasta
done

for a in *.RLP.LysM.fasta
do
cat $a >>onion.RLP.LysM.fasta
done

for a in *.TMCC.all.fasta
do
cat $a >>onion.TMCC.all.fasta
done

for a in *_slrk.fasta
do
cat $a >>onion_slrk.fasta
done




