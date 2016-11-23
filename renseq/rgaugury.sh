#!/bin/bash
input=/home/sobczm/popgen/renseq/input/transcriptomes
scripts=/home/sobczm/bin/popgen/renseq
#Additional tool to annotate non-canonical R genes not containing the NBS domain.
#Using it to annotate the following groups of candidate genes:
#RLK (transmembrane-LRR/LysM-STTK domain),
#RLP (transmembrane-LRR/LysM domain)
#TMCC (transmembrane-coiled coil)

#Create a seperate directory structure
cd $input
mkdir -p rgaugury/kim
cp ./KIM/onion_kim_protein.fa ./rgaugury/kim
mkdir -p rgaugury/raj
cp ./RAJ/onion_raj_protein.fa ./rgaugury/raj
mkdir -p rgaugury/nz
cp ./NZ/onion_nz_protein.fa ./rgaugury/nz
mkdir -p rgaugury/maria
cp ./MARIA/onion_maria_protein.fa ./rgaugury/maria
mkdir -p rgaugury/brian
cp ./BRIAN/onion_brian_protein.fa ./rgaugury/brian

#Remove stop-codons from six-frame translated contigs
sed -i 's/\*//g' ./rgaugury/kim/onion_kim_protein.fa
sed -i 's/\*//g' ./rgaugury/raj/onion_raj_protein.fa
sed -i 's/\*//g' ./rgaugury/nz/onion_nz_protein.fa
sed -i 's/\*//g' ./rgaugury/brian/onion_brian_protein.fa
sed -i 's/\*//g' ./rgaugury/maria/onion_maria_protein.fa

#Break down the files into smaller chunks and run rgaugury for each transcriptome.
#As this taking too long, going to parallelise in a crude way by splitting
#into small files with 100 protein sequences and running 20 at a time.
cd $input/rgaugury/kim/
file=onion_kim_protein.fa

#!/bin/bash
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $file>temp && mv temp $file
awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%100==0){file=sprintf("myseq%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < $file
for file in myseq*.fa
do
    Jobs=$(qstat | grep 'sub_rgaugu' | wc -l)
    while [ $Jobs -gt 40 ]
    do
        sleep 10
        printf "."
        Jobs=$(qstat | grep 'sub_rgaugu' | wc -l)
    done
qsub $scripts/sub_rgaugury.sh $file
done

#repeat the above mini script but changing the dirs and input files for other transcriptomes
cd $input/rgaugury/raj/
file=onion_raj_protein.fa
cd $input/rgaugury/nz/
file=onion_nz_protein.fa
cd $input/rgaugury/brian/
file=onion_brian_protein.fa
cd $input/rgaugury/maria/
file=onion_maria_protein.fa
