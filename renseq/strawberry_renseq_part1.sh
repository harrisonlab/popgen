#!/bin/bash
scripts=/home/sobczm/bin/popgen/renseq
input=/home/sobczm/popgen/renseq/strawberry

cd $input
mkdir genome
cd genome

#Ver 1.1 
cp /home/vicker/octoploid_mapping/vesca/vesca_assemblies/v1.1/* ./
#The 2015 Towson Uni annotation
wget ftp://ftp.bioinfo.wsu.edu/species/Fragaria_vesca/Fvesca-genome.v1.1.a2/genes/Fragaria_vesca_v1.1.a2.gff3.gz
#Combine unanchored contigs with all chromosomes
cat fvesca_v1.1_pseudo.fna >fvesca_v1.1_all.fa
cat fvesca_v1.1_unanchored.fna >>fvesca_v1.1_all.fa
#Need to correct the GFF file for the unanchored chromosome location to be given as LG0
sed -i 's/LGUn/LG0/g' Fragaria_vesca_v1.1.a2.gff3 
gffread=/home/sobczm/bin/gffread/gffread/gffread 
$gffread -w fvesca_v1.1_all_annotated.fa -g fvesca_v1.1_all.fa Fragaria_vesca_v1.1.a2.gff3

#Translate the CDS sequences
emboss=/home/armita/prog/emboss/EMBOSS-4.0.0/bin
$emboss/transeq fvesca_v1.1_all_annotated.fa -outseq fvesca_v1.1_all_proteins.fa

#Remove stop codons from protein ends
a=fvesca_v1.1_all_proteins.fa
b=$( echo $a | sed -e 's/.fa/\_nostop.fa/' )
sed 's/\*//g' $a > $b

#Run Rgaugury on pep sequences 
mkdir -p $input/rgaugury/vesca1.1
cp $input/genome/fvesca_v1.1_all_proteins_nostop.fa $input/rgaugury/vesca1.1

cd $input/rgaugury/vesca1.1

file=fvesca_v1.1_all_proteins_nostop.fa
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $file>temp && mv temp $file
awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%100==0){file=sprintf("myseq%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < $file

for file in myseq*.fa
do
    Jobs=$(qstat | grep 'sub_rgaugu' | wc -l)
    while [ $Jobs -gt 100 ]
    do
        sleep 10
        printf "."
        Jobs=$(qstat | grep 'sub_rgaugu' | wc -l)
    done
qsub $scripts/sub_rgaugury.sh $file
done
done
