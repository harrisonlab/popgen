#!/bin/bash
scripts=/home/sobczm/bin/popgen/renseq
input=/home/sobczm/popgen/renseq/apple_2017

#Re-prediction of R genes and new bait design using the latest Golden Delicious genome from Daccord et al. (2017)

#Download genome, genome annotation and proteins.
wget https://iris.angers.inra.fr/gddh13/downloads/GDDH13_1-1_formatted.fasta.bz2
wget https://iris.angers.inra.fr/gddh13/downloads/gene_models_20170612.gff3.bz2
wget https://iris.angers.inra.fr/gddh13/downloads/GDDH13_1-1_prot.fasta.bz2

#Extract CDS from mRNA.
awk '$3 == "CDS" {print $0}' gene_models_20170612.gff3 >gene_models_20170612_cds.gff3
##Extract the sequences of the genes to be used as blast database (CDS + 5'/3' UTR)
genome_seq=GDDH13_1-1_formatted.fasta
gffread=/home/sobczm/bin/gffread/gffread/gffread
$gffread -w GDDH13_1-1_formatted_cds.fasta -g $genome_seq gene_models_20170612_cds.gff3
sed -i 's/mRNA://g' GDDH13_1-1_formatted_cds.fasta 
#Split the protein sequences into small files with 100 protein sequences and run 50 at a time.
file=GDDH13_1-1_prot.fasta
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

#The plan is to have 1 bait library of NBS-LRR genes (20,000 baits) and one of RLK and RLP genes (20,000 baits).