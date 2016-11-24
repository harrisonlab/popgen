#!/bin/bash
input=/home/sobczm/popgen/renseq/input/transcriptomes
input2=/home/sobczm/popgen/renseq/input/transcriptomes/plant_rgenes/common_onion
scripts=/home/sobczm/bin/popgen/renseq
rgenes=/home/sobczm/bin/plant_rgenes

#Download the following files from the Pfam ftp site (Pfam v. 30.0)
#Pfam-A.hmm
#Pfam-A.hmm.dat
#active_site.dat

#Generate binary files for Pfam-A.hmm by running the following commands:
hmmpress /home/sobczm/bin/hmmer-3.1b2/pfam30/Pfam-A.hmm

#Prepare 6 frame translations of all contigs in the assemblies
java -jar $scripts/Translate6Frame.jar -i $input/KIM/GBRQ01_1_fsa_nt_combined_kim.fasta \
-o $input/KIM/onion_kim_protein.fa
java -jar $scripts/Translate6Frame.jar -i $input/NZ/GBGJ01_1_fsa_nt_nz.fasta \
-o $input/NZ/onion_nz_protein.fa
java -jar $scripts/Translate6Frame.jar -i $input/RAJ/GBJZ01_1_fsa_nt_raj.fasta \
-o $input/RAJ/onion_raj_protein.fa
java -jar $scripts/Translate6Frame.jar -i $input/NZ/GBGJ01_1_fsa_nt_nz.fasta \
-o $input/NZ/onion_nz_protein.fa
java -jar $scripts/Translate6Frame.jar -i $input/RAJ/GBJZ01_1_fsa_nt_raj.fasta \
-o $input/RAJ/onion_raj_protein.fa
java -jar $scripts/Translate6Frame.jar -i $input/MARIA/Trinity.fasta \
-o $input/MARIA/onion_maria_protein.fa
java -jar $scripts/Translate6Frame.jar -i $input/BRIAN/Trinity_11Oct2016.fasta \
-o $input/BRIAN/onion_brian_protein.fa

#Execute the pipeline using all 6-frame translations of the assemblies
#Substitute the run_pfam_scan.sh with a simplified command
pfam_scan.pl -e_seq 1 -e_dom 1 -as -outfile onion_kim_protein.out -cpu 8 -fasta $input/KIM/onion_kim_protein.fa \
-dir /home/sobczm/bin/hmmer-3.1b2/pfam30

#As this taking too long, going to parallelise in a crude way by splitting
#into small files with 100 protein sequences and running 20 at a time.

names=( "kim" "raj" "nz" "brian" "maria" )

for f in "${names[@]}"
do
fca=$(echo $f | tr '[:lower:]' '[:upper:]')
file=onion_${f}_protein.fa
mkdir -p $input2/$fca
cp $input/$fca/$file $input2/$fca/$file
cd $input2/$fca
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $file>temp && mv temp $file
awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%100==0){file=sprintf("myseq%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < $file
for file in myseq*.fa
do
    Jobs=$(qstat | grep 'sub_hmmsca' | wc -l)
    while [ $Jobs -gt 100 ]
    do
        sleep 10
        printf "."
        Jobs=$(qstat | grep 'sub_hmmsca' | wc -l)
    done
qsub $scripts/sub_hmmscan.sh $file
done
#concatenate the outputs into one file
for a in $input2/$fca/*.out; do cat $a >> onion_${f}_protein.out; done
#Prepare input for domain analysis
mkdir -p $input2/domains/$f
cp $input2/$fca/onion_${f}_protein* $input2/domains/$f
cp $input2/test2/db_descriptions.txt $input2/domains/$f
mv $input2/domains/$f/onion_${f}_protein.fa onion_167_TAIR10.protein.fa
mv $input2/domains/$f/onion_${f}_protein.out onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.out
#Carry out domain parsing
perl $rgenes/processing_scripts/K-parse_Pfam_domains_v3.1.pl --pfam $input2/domains/$f/onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.out \
--evalue 0.0001 --out $input2/domains/$f/onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.parsed.verbose --verbose T
perl $rgenes/processing_scripts/K-parse_Pfam_domains_NLR-fusions-v2.3.pl \
--indir $input2/domains/$f --evalue 0.0001 --outdir $input2/domains/$f--db_description $input/domains/$f/db_descriptions.txt
done

#Sort the output
for a in $input/domains/*/onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.parsed.verbose.NLR*
do
sort -k 1 $a >${a%.*}_sorted.txt
done
