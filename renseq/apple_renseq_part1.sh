#!/bin/bash
scripts=/home/sobczm/bin/popgen/renseq
input=/home/sobczm/popgen/renseq/apple

#Download all the apple genome resources
cd $input
mkdir genome
cd genome
wget -r ftp://climb.genomics.cn/pub/10.5524/100001_101000/100189/
#File with CDS sequences
ls -l $input/genome/EVM.out.cds
#File with pep sequences
ls -l $input/genome/gene/EVM.out.pep

#Download the ver 3.0.a1 of the original apple genome 
#CDS and pep
cd $input/genome
mkdir Velasco && cd Velasco
wget ftp://ftp.bioinfo.wsu.edu/species/Malus_x_domestica/Malus_x_domestica-genome.v3.0.a1/genes/Malus_x_domestica.v3.0.a1_gene_set_cds.fasta.gz
wget ftp://ftp.bioinfo.wsu.edu/species/Malus_x_domestica/Malus_x_domestica-genome.v3.0.a1/genes/Malus_x_domestica.v3.0.a1_gene_set_pep.fasta.gz
wget ftp://ftp.bioinfo.wsu.edu/species/Malus_x_domestica/Malus_x_domestica-genome.v3.0.a1/assembly/Malus_x_domestica.v3.0.a1_contigs.gff3.gz
wget ftp://ftp.bioinfo.wsu.edu/species/Malus_x_domestica/Malus_x_domestica-genome.v3.0.a1/assembly/Malus_x_domestica.v3.0.a1_contigs.fasta.gz
wget ftp://ftp.bioinfo.wsu.edu/species/Malus_x_domestica/Malus_x_domestica-genome.v3.0.a1/genes/Malus_x_domestica.v3.0.a1_v1_gene_alignemnt.gff3.gz

#Remove terminal codons (*) from the pep file, as InterProScan
#does not accept sequences containing stop codons.
a=Malus_x_domestica.v3.0.a1_gene_set_pep.fasta
b=$( echo $a | sed -e 's/.fasta/\_nostop.fasta/' )
sed 's/\*//g' $a > $b
#Run Rgaugury on pep sequences in both genomes
mkdir -p $input/rgaugury/Velasco
cp $input/genome/Velasco/Malus_x_domestica.v3.0.a1_gene_set_pep_nostop.fasta $input/rgaugury/Velasco
mkdir -p $input/rgaugury/Li
cp $input/genome/gene/EVM.out.pep $input/rgaugury/Li
#Split the sequences into small files with 100 protein sequences and run 50 at a time.
cd $input/rgaugury/Velasco
file=Malus_x_domestica.v3.0.a1_gene_set_pep_nostop.fasta
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

#Run Phobius and final prediction
rm *RLKorRLP.domain.prediction.txt
rgaugury=/home/sobczm/bin/rgaugury/RGAugury.pl
for file in myseq*.fa
do
perl $rgaugury -p $file 
done

cd $input/rgaugury/Velasco
mkdir final
mv myseq* ./final

cd $input/rgaugury/Li 
file=EVM.out.pep
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

#Run Phobius and final prediction
rm *RLKorRLP.domain.prediction.txt
rgaugury=/home/sobczm/bin/rgaugury/RGAugury.pl
for file in myseq*.fa
do
perl $rgaugury -p $file 
done

cd $input/rgaugury/Li 
mkdir final
mv myseq* ./final