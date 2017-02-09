#!/bin/bash
scripts=/home/sobczm/bin/popgen/renseq
input=/home/sobczm/popgen/renseq/apple

#Extract select contigs in each assembly

#Li
cd $input/rgaugury/Li/lists
for list in *.no_transposon
do
python $scripts/keep_list_genes2.py $list $input/genome/gene/EVM.out.cds No
done
cp *.fasta $input/baits

#Velasco
cd $input/rgaugury/Velasco/lists
for list in *.no_transposon
do
python $scripts/keep_list_genes2.py $list $input/genome/Velasco/Malus_x_domestica.v3.0.a1_gene_set_cds.fasta No
done
cp *.fasta $input/baits

#Concatenate the contig sequence files for each gene family
cd $input/baits

for a in *.NBS.all.fasta
do
cat $a >>apple.NBS.all.fasta
done

for a in *.NBS.CN.fasta
do
cat $a >>apple.NBS.CN.fasta
done

for a in *.NBS.CNL.fasta
do
cat $a >>apple.NBS.CNL.fasta
done

for a in *.NBS.NBS.fasta
do
cat $a >>apple.NBS.NBS.fasta
done

for a in *.NBS.NL.fasta
do
cat $a >>apple.NBS.NL.fasta
done

for a in *.NBS.OTHER.fasta
do
cat $a >>apple.NBS.OTHER.fasta
done

for a in *.NBS.TN.fasta
do
cat $a >>apple.NBS.TN.fasta
done

for a in *.NBS.TNL.fasta
do
cat $a >>apple.NBS.TNL.fasta
done

for a in *.NBS.TX.fasta
do
cat $a >>apple.NBS.TX.fasta
done

for a in *.RLK.all.fasta
do
cat $a >>apple.RLK.all.fasta
done

for a in *.RLK.LRR.fasta
do
cat $a >>apple.RLK.LRR.fasta
done

for a in *.RLK.LysM.fasta
do
cat $a >>apple.RLK.LysM.fasta
done

for a in *.RLK.other.fasta
do
cat $a >>apple.RLK.other.fasta
done

for a in *.RLP.all.fasta
do
cat $a >>apple.RLP.all.fasta
done

for a in *.RLP.LRR.fasta
do
cat $a >>apple.RLP.LRR.fasta
done

for a in *.RLP.LysM.fasta
do
cat $a >>apple.RLP.LysM.fasta
done

for a in *.TMCC.all.fasta
do
cat $a >>apple.TMCC.all.fasta
done

for a in *_slrk.fasta
do
cat $a >>apple_slrk.fasta
done