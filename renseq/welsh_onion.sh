#!/bin/bash
input=/home/sobczm/popgen/renseq/input/reads
scripts=/home/sobczm/bin/popgen/renseq

#As the number of R genes discovered in the onion transcriptomes was quite low
#mining the other Allium species transcriptomes for them. Towards this end
#will download the reads for published Welsh onion (and also remaining common onion)
#transcriptomes and re-assamble them with Trinity.
cd $input
sh ./download_reads.sh

for reads in $input/*/*/all_reads_trim_1.fastq.gz
do
    Jobs=$(qstat | grep 'sub_trinit' | wc -l)
    while [ $Jobs -gt 1 ]
    do
        output_dir=$(echo $reads | awk -F/ '{print $(NF-1)}')
        echo $output_dir
        sleep 100
        printf "."
        Jobs=$(qstat | grep 'sub_trinit' | wc -l)
    done
reads2=$(echo "$reads" | sed 's/_1.fastq.gz/_2.fastq.gz/')
qsub $scripts/sub_trinity_assembly.sh $reads $reads2 $output_dir
done

#Following that, carry on with the previously identified analyses to uncover
#R genes in the transcriptomes: NLR Parser, Plant R Genes and rgaugury
input=/home/sobczm/popgen/renseq/input/transcriptomes
input2=/home/sobczm/popgen/renseq/input/transcriptomes/plant_rgenes/common_onion
scripts=/home/sobczm/bin/popgen/renseq

##Common onion
#1) NLRParser
cd $input/CORNELL
sh $scripts/sub_nlrparser.sh cornell_Trinity.fasta
sort -k 1 cornell_Trinity_nlr.tsv >cornell_Trinity_nlr_sorted.tsv
cd $input/H6
sh $scripts/sub_nlrparser.sh h6_Trinity.fasta
sort -k 1 h6_Trinity_nlr.tsv >h6_Trinity_nlr_sorted.tsv
cd $input/HAN
sh $scripts/sub_nlrparser.sh han_Trinity.fasta
sort -k 1 han_Trinity_nlr.tsv >han_Trinity_nlr_sorted.tsv
cd $input/SP3B
sh $scripts/sub_nlrparser.sh sp3b_Trinity.fasta
sort -k 1 sp3b_Trinity_nlr.tsv >sp3b_Trinity_nlr_sorted.tsv
cd $input/SP3B
sh $scripts/sub_nlrparser.sh sp3b_Trinity.fasta
sort -k 1 sp3b_Trinity_nlr.tsv >sp3b_Trinity_nlr_sorted.tsv
cd $input/LIU
sh $scripts/sub_nlrparser.sh liu_Trinity.fasta
sort -k 1 liu_Trinity_nlr.tsv >liu_Trinity_nlr_sorted.tsv
cd $input/SUN
sh $scripts/sub_nlrparser.sh sun_Trinity.fasta
sort -k 1 sun_Trinity_nlr.tsv >sun_Trinity_nlr_sorted.tsv

#Generate 6 frame protein translations
java -jar $scripts/Translate6Frame.jar -i $input/CORNELL/cornell_Trinity.fasta \
-o $input/CORNELL/onion_cornell_protein.fa
java -jar $scripts/Translate6Frame.jar -i $input/H6/h6_Trinity.fasta \
-o $input/H6/onion_h6_protein.fa
java -jar $scripts/Translate6Frame.jar -i $input/HAN/han_Trinity.fasta \
-o $input/HAN/onion_han_protein.fa
java -jar $scripts/Translate6Frame.jar -i $input/SP3B/sp3b_Trinity.fasta \
-o $input/SP3B/onion_sp3b_protein.fa
java -jar $scripts/Translate6Frame.jar -i $input/LIU/liu_Trinity.fasta \
-o $input/LIU/onion_liu_protein.fa
java -jar $scripts/Translate6Frame.jar -i $input/SUN/sun_Trinity.fasta \
-o $input/SUN/onion_sun_protein.fa

#2) Plant R genes pipeline
names=( "cornell" "h6" "han" "sp3b" "liu" "sun" )
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
