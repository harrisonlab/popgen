#!/bin/bash
wdir=/home/sobczm/popgen/codon/blast
scripts=/home/sobczm/bin/popgen/codon
#2) filter blast results by percentage id and subject coverage
cd $wdir
for a in *one*.db
do
python $scripts/filter_blast.py --i $a
done

for z in *filtered
do
cp $z ./dagchainer
done

cd ./dagchainer
#3) blast_to_dagchainer.py
